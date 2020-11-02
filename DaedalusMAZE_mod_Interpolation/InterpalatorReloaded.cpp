/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <iostream>
#include "boost/multi_array.hpp"
#include "common.h"
#include "io.h"
#include "interpolation.h"
#include "cuInterp.h"
#include <fstream>
#include <future>


int main(int argc, char *argv[]){

    time_t start = time(0);
    int ncerror;
    static const int NC_ERR = 2;

    /* Define and Start Allocating Memory*/
    std::vector<float> nextVar, nextAlts,nextGtime;  //used for preloading files
    std::vector<int> filesToRead;
    std::vector<float> glat,glon,glev,gtime;
    std::vector<float> dalt,dlon,dlat,acc1,acc2,acc3,intdata,dlatRun, dlonRun, daltRun, intdataRun,dtime,var,alts;
    bool nextFileLoaded = false;
    std::string varName;
    bool ferror;
   
    /*Write Init Log*/
    logInit(argv,argc);
    std::ofstream logger;
    logger.open("logfile.txt", std::ios_base::app);

    /*Read Orbit And Domain Grid*/
    ncerror = ReadGrid(glat,glon,glev,argv[1]);
    ferror=buildDtime(dtime);
    ncerror = ReadOrbit(dlon,dlat,dalt,acc1,acc2,acc3);
    ferror=WriteOrbit(dlat,dlon,dalt,dtime,acc1,acc2,acc3);


    /*Split Orbit into proper segments*/
    const int numPasses = orbit::DIM / orbit::BATCH;
    int leftoverOrbit =orbit::DIM%orbit::BATCH;
    int leftoverFiles =(argc-1)%model::BATCH;
    //std::vector<std::string> vars={"ZGMID"};    


    std::vector<std::string> vars = {"BX_si", "BY_si", "BZ_si", "EEX_si", "EEY_si", "EEZ_si", "SIGMA_PED","SIGMA_HAL", "UN_si", "VN_si", "Wi_lev_si", "Ohmic","ZGMID","Wind_heating","Convection_heating","DEN","QJOULE"};
    /*Sanity Check ()Issue #2)*/
    if (numPasses > argc){
      std::cout<<"Aborting at Line:"<< __LINE__  << "\tOrbit Size probably does not match model data file count [argc]\t"  <<std::endl; abort();
    }

    time_t interval0 = time(0);

    for (int numVar=0; numVar<vars.size();++numVar){  /*variable loop*/
        varName=vars.at(numVar);
        

        for (int i=0; i<numPasses;++i){ /*file loop*/


            logger << "Pass:"<<i<<"\t"<< varName<<std::endl;


            /*Read Batch*/
            if (RUN::MANUAL){
                filesToRead.push_back(1);
            }else{
                for (int numberOfFile=0; numberOfFile<model::BATCH; ++numberOfFile){
                    filesToRead.push_back(i*model::BATCH+1+numberOfFile);
                    logger <<"\t"<< argv[i*model::BATCH+1+numberOfFile] << std::endl;
                }
            }


            /*Read Variables*/
            if (!nextFileLoaded){
                ncerror = ReadVar(var, alts, gtime, argv, varName, filesToRead);
                if (ncerror == NC_ERR){  std::cout<<"Aborting at Line:"<< __LINE__  <<std::endl; abort();}
            }
            
            #ifdef PRELOAD_FILES
                #pragma message( "PRELOADING FILES IS ENABLED.USE WITH CAUTION ")
                // Check if PreLoading is Possible
                std::array<int,3> currentState= {numVar,i,numPasses}; /* current state to estimate if preloading is available */
                bool preLoadFiles = PreaLoadingAvailable(currentState);
                std::future<bool> fut = std::async(std::launch::async,ReadAsyncVar,std::ref(nextVar), std::ref(nextAlts), std::ref(nextGtime),argv,varName, filesToRead.at(0)+1,preLoadFiles);   
            #endif

            /*Get Orbit Batch*/
            for (int j=0; j<orbit::BATCH;++j){
                dlatRun.push_back(dlat.at(i*orbit::BATCH+j));
                dlonRun.push_back(dlon.at(i*orbit::BATCH+j));
                daltRun.push_back(dalt.at(i*orbit::BATCH+j));                
            }

            /*Interpolate*/
            if (RUN::DEVICE=="CPU"){  /*CPU interpolation*/

                if (RUN::METHOD=="TRILINEAR"){    
                    ferror=InterpolateBatch(glat,glon,glev,gtime,dlatRun,dlonRun,daltRun,alts,var,intdataRun);
                }
                if (RUN::METHOD=="TRICUBIC"){
                    ferror = TriCubicSpline(glat, glon, glev, gtime, dlatRun, dlonRun, daltRun, alts, var, intdataRun,vars.at(numVar));
                //    AdaptiveLPF(intdataRun,daltRun,1);
                }
                
            }else{  /*GPU Device interpolation*/
                #ifdef GPU
                ferror=PrepeareForGPU(glat,glon,glev,gtime,dlatRun,dlonRun,daltRun,alts,var,intdataRun);
                #endif
            }
            
            
            /*Store Batch*/
            intdata.insert(intdata.end(), intdataRun.begin(), intdataRun.end());

            
            /*Reset*/
            dlatRun.resize(0);
            dlonRun.resize(0);
            daltRun.resize(0);
            intdataRun.resize(0);
            filesToRead.resize(0);

            #ifdef PRELOAD_FILES
                // Grab Preloaded File
                nextFileLoaded = fut.get();
                if (nextFileLoaded){
                    alts=nextAlts;
                    var=nextVar;
                    gtime=nextGtime;
                }
            #endif
            
        }//End Of Main Loop


        /* Intrpolate leftOver files And Orbit*/
        if (leftoverFiles>0 && leftoverOrbit>0){
            logger << "\t" << "LeftOver Orbit"<< std::endl;
            for (int kk = dlat.size()-leftoverOrbit; kk < orbit::DIM; ++kk){
                dlatRun.push_back(dlat.at(kk));
                dlonRun.push_back(dlon.at(kk));
                daltRun.push_back(dalt.at(kk));
            }

            if (RUN::MANUAL){
                filesToRead.push_back(1);
            }else{

                for (int kk = argc-leftoverFiles; kk < argc; ++kk){
                    filesToRead.push_back(kk);
                    logger << "\t" << filesToRead.at(kk)<< std::endl;
                }
            }
            

            if (RUN::METHOD=="TRILINEAR"){
                ferror = InterpolateBatch(glat, glon, glev, gtime, dlatRun, dlonRun, daltRun, alts, var, intdataRun);
            }
            if (RUN::METHOD=="TRICUBIC"){
                ferror = TriCubicSpline(glat, glon, glev, gtime, dlatRun, dlonRun, daltRun, alts, var, intdataRun,vars.at(numVar) );
            }
            

            //Merge Data
            intdata.insert(intdata.end(), intdataRun.begin(), intdataRun.end());

            /*Reset*/
            dlatRun.resize(0);
            dlonRun.resize(0);
            daltRun.resize(0);
            intdataRun.resize(0);
            filesToRead.resize(0);

        }
      
    
        WriteData(intdata, varName);  
        intdata.resize(0);
        logger << "Saving Data" << std::endl;
        time_t interval1 = time(0);
        double predictedTime = difftime(interval1, interval0);
        logger<< "ETA=\t"<< predictedTime*(vars.size()-numVar)/60.<<"\t minutes"<<std::endl;
    }
    
 ;

    time_t end = time(0);
    double time = difftime(end, start);
    printf("Done in %f minutes\n", time/60.0);
    logger << "Interpolator Done in :"<<time/60.0<<"minutes" << std::endl;
    logger.close();
}
