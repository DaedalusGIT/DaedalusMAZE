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

#include <iostream>
#include <netcdf>
#include <stdio.h>
#include "common.h"
#include <omp.h>
#include "interpolation.h"
#include<array>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



bool ReadGrid(std::vector<float> &glatV,
            std::vector<float> &glonV,
            std::vector<float> &glevV ,
            char *sampleFile          

)

{
    auto glat = new float[grid::NLAT]();
    auto glon = new float[grid::NLON]();
    auto glev = new float[grid::NALT]();

    static const int NC_ERR = 2;

    NcFile dataFile(sampleFile, NcFile::read);
    NcVar latVar, lonVar,levVar,timeVar;
    latVar = dataFile.getVar("lat");
     if(latVar.isNull()) return NC_ERR;
    lonVar = dataFile.getVar("lon");
     if(lonVar.isNull()) return NC_ERR;
    levVar = dataFile.getVar("lev");
    if (levVar.isNull())
        return NC_ERR;

    latVar.getVar(glat);
    lonVar.getVar(glon);
    levVar.getVar(glev);

    for (int i = 0; i<grid::NLAT; ++i){
        glatV.push_back(glat[i]);
    }

    for (int i = 0; i<grid::NLON; ++i){
        glonV.push_back(glon[i]);
    }

    for (int i = 0; i<grid::NALT; ++i){
        glevV.push_back(glev[i]);
    }

    delete[] glat;
    delete[] glon;
    delete[] glev;

    return true;
}

bool ReadVar(std::vector<float> &varVector,
            std::vector<float> &altVector,
            std::vector<float> &gtimeV,
            char *files[],
            std::string varName,
            std::vector<int> numFile)
{
    static const int NC_ERR = 2;
    auto var = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
    auto alts = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
    auto gtime = new float[grid::NTIMES]();
    std::string filename;
    varVector.resize(0);
    altVector.resize(0);
    gtimeV.resize(0);


    for (int start=0; start<numFile.size();++start){

        filename=files[numFile.at(start)];

        NcFile dataFile(filename, NcFile::read);
        NcVar timeVar,altVar,dataVar;
        dataVar = dataFile.getVar(varName);
        if (dataVar.isNull()) return NC_ERR;
        altVar = dataFile.getVar("ZGMID");
        if (altVar.isNull()) return NC_ERR;
        timeVar = dataFile.getVar("time");
        if (timeVar.isNull()) return NC_ERR;
        altVar.getVar(alts);
        dataVar.getVar(var);
        timeVar.getVar(gtime);
        

        #pragma omp collapse(4)    
        for (int z=0; z<grid::NLON; ++z){
            for (int k=0; k<grid::NLAT; ++k){
                for (int j=0; j<grid::NALT; ++j){
                    for (int i=0; i<grid::NTIMES; ++i){

                        varVector.push_back(var[i][j][k][z]);
                        altVector.push_back(alts[i][j][k][z]);
        
        
                    }
                }
            }
        }

        for (int i=0; i<grid::NTIMES; ++i){
            gtimeV.push_back(gtime[i]);
        }   

    }

    delete [] var;
    delete [] alts;
    delete [] gtime;
    // std::cout<<gtimeV.size()<<std::endl;

    return true;
}


bool ReadAsyncVar(std::vector<float> &varVector,
            std::vector<float> &altVector,
            std::vector<float> &gtimeV,
            char *files[],
            std::string varName,
            int FileID,
            bool filesExist)
{

    if (!filesExist){
        return false;
    }
    static const int NC_ERR = 2;
    auto var = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
    auto alts = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
    auto gtime = new float[grid::NTIMES]();
    std::string filename;
    varVector.resize(0);
    altVector.resize(0);
    gtimeV.resize(0);


    

    filename=files[FileID];
    NcFile dataFile(filename, NcFile::read);
    NcVar timeVar,altVar,dataVar;
    dataVar = dataFile.getVar(varName);
    if (dataVar.isNull()) return NC_ERR;
    altVar = dataFile.getVar("ZGMID");
    if (altVar.isNull()) return NC_ERR;
    timeVar = dataFile.getVar("time");
    if (timeVar.isNull()) return NC_ERR;
    altVar.getVar(alts);
    dataVar.getVar(var);
    timeVar.getVar(gtime);
    

    #pragma omp collapse(4)    
    for (int z=0; z<grid::NLON; ++z){
        for (int k=0; k<grid::NLAT; ++k){
            for (int j=0; j<grid::NALT; ++j){
                for (int i=0; i<grid::NTIMES; ++i){

                    varVector.push_back(var[i][j][k][z]);
                    altVector.push_back(alts[i][j][k][z]);
    
    
                }
            }
        }
    }

    for (int i=0; i<grid::NTIMES; ++i){
        gtimeV.push_back(gtime[i]);
    }   



    delete [] var;
    delete [] alts;
    delete[] gtime;
    // std::cout<<gtimeV.size()<<std::endl;

    return true;
}


bool ReadMultiVar(std::vector<float> &varVector,
             std::vector<float> &altVector,
             std::vector<float> &gtimeV,
             char *files[],
             std::string varName,
             std::vector<int> numFile)
 {
     static const int NC_ERR = 2;
     auto bigvar = new float[grid::NTIMES*model::BATCH][grid::NALT][grid::NLAT][grid::NLON]();
     auto bigalts = new float[grid::NTIMES*model::BATCH][grid::NALT][grid::NLAT][grid::NLON]();
     auto biggtime = new float[grid::NTIMES*model::BATCH]();

     auto var = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
     auto alts = new float[grid::NTIMES][grid::NALT][grid::NLAT][grid::NLON]();
     auto gtime = new float[grid::NTIMES]();
     std::string filename;

     varVector.resize(0);
     altVector.resize(0);
     gtimeV.resize(0);

     for (int start=0; start<numFile.size();++start){


        
         filename=files[numFile.at(start)];

         NcFile dataFile(filename, NcFile::read);
         NcVar timeVar,altVar,dataVar;
         dataVar = dataFile.getVar(varName);
         if (dataVar.isNull()) return NC_ERR;
         altVar = dataFile.getVar("ZGMID");
         if (altVar.isNull()) return NC_ERR;
         timeVar = dataFile.getVar("time");
         if (timeVar.isNull()) return NC_ERR;
       
         altVar.getVar(alts);
         dataVar.getVar(var);
         timeVar.getVar(gtime);


         int pos1=start*grid::NTIMES;
         int pos2=start*grid::NTIMES+grid::NTIMES;

         int smalltime = 0;
         for (int i=pos1; i<pos2; ++i){
             for (int j=0; j<grid::NALT; ++i){
                 for (int k=0; k<grid::NLAT; ++i){
                     for (int z=0; z<grid::NLON; ++i){

                        bigvar[i][j][k][z]=var[smalltime][j][k][z];
                        bigalts[i][j][k][z]=alts[smalltime][j][k][z];
                   
                   
                     }
                 }
             }
         smalltime+=1;
         }


         for (int i=0; i<grid::NTIMES; ++i){
             gtimeV.push_back(gtime[i]);
         }
                    
     }




         #pragma omp collapse(4)    
         for (int z=0; z<grid::NLON; ++z){
             for (int k=0; k<grid::NLAT; ++k){
                 for (int j=0; j<grid::NALT; ++j){
                     for (int i=0; i<grid::NTIMES*model::BATCH; ++i){

                         varVector.push_back(bigvar[i][j][k][z]);
                         altVector.push_back(bigalts[i][j][k][z]);
        
        
                     }
                 }
             }
         }


     delete [] var;
     delete [] alts;
     delete[] gtime;
     delete [] bigvar;
     delete [] bigalts;
     delete[] biggtime;
     // std::cout<<gtimeV.size()<<std::endl;

     return true;
 }


bool ReadOrbit(std::vector<float> &lonV,
               std::vector<float> &latV,
               std::vector<float> &altV,
               std::vector<float> &acc1,
               std::vector<float> &acc2,
               std::vector<float> &acc3)
{
    static const int NC_ERR = 2;
    auto lon = new float [orbit::DIM]();
    auto lat = new float [orbit::DIM]();
    auto alt = new float [orbit::DIM]();
    auto mlt1 = new float [orbit::DIM]();
    auto mlt2 = new float [orbit::DIM]();
    auto mlt3 = new float [orbit::DIM]();
    
    std::string name ="lat";
    NcFile dataFile(orbit::OrbitName, NcFile::read);
    // NcFile dataFile("../Input/2.nc", NcFile::read);
    NcVar lonVar, altVar, latVar,ext1,ext2,ext3;
    latVar = dataFile.getVar(name);
    if (latVar.isNull()) return NC_ERR;
    altVar = dataFile.getVar("altitude");
    if (altVar.isNull())return NC_ERR;
    lonVar = dataFile.getVar("lon");
    if (lonVar.isNull()) return NC_ERR;
    ext1 = dataFile.getVar("DaedalusMLT");
    if (ext1.isNull()) return NC_ERR;
    ext2 = dataFile.getVar("DaedalusMagneticLongitude");
    if (ext2.isNull()) return NC_ERR;
    ext3 = dataFile.getVar("DaedalusMagneticLatitude");
    if (ext3.isNull()) return NC_ERR;
    altVar.getVar(alt);
    latVar.getVar(lat);
    lonVar.getVar(lon);
    ext1.getVar(mlt1);
    ext2.getVar(mlt2);
    ext3.getVar(mlt3);

    for (int i=0; i<orbit::DIM;++i){
        latV.push_back(lat[i]);
        if (lon[i]>=180.0){
            lonV.push_back(lon[i] - 360.0);
        }
        else{
            lonV.push_back(lon[i]);
        }
        altV.push_back(alt[i]);
        acc1.push_back(mlt1[i]);
        acc2.push_back(mlt2[i]);
        acc3.push_back(mlt3[i]);

    }


    delete[] lon;
    delete[] alt;
    delete[] lat;
    delete[] mlt1;
    delete[] mlt2;
    delete[] mlt3;

    return true;
}

bool WriteOrbit(std::vector<float> dlat,
                std::vector<float> dlon,
                std::vector<float> dalt,
                std::vector<float> dtime,
                std::vector<float> acc1,
                std::vector<float> acc2,
                std::vector<float> acc3)
{

    NcFile dataFile(orbit::OutFile, NcFile::replace);

    auto outtime = new float[dlat.size()];
    auto outlat = new float[dlat.size()];
    auto outlon = new float[dlat.size()];
    auto outalt = new float[dlat.size()];
    auto outacc1 = new float[dlat.size()];
    auto outacc2 = new float[dlat.size()];
    auto outacc3 = new float[dlat.size()];
    
   
    
    for (int i = 0; i < dlat.size(); ++i)
    {

        outtime[i] = dtime.at(i);
        outlat[i] = dlat.at(i);
        outlon[i] = dlon.at(i);
        outalt[i] = dalt.at(i);
        outacc1[i] = acc1.at(i);
        outacc2[i] = acc2.at(i);
        outacc3[i] = acc3.at(i);
    }

    vector<NcDim> dims;
    NcDim timeDim = dataFile.addDim("time", dtime.size());
    dims.push_back(timeDim);

    NcVar lat = dataFile.addVar("lat", ncFloat, dims);
    NcVar time = dataFile.addVar("time", ncFloat, dims);
    NcVar lon = dataFile.addVar("lon", ncFloat, dims);
    NcVar alt = dataFile.addVar("alt", ncFloat, dims);
    NcVar mlt1 = dataFile.addVar("DaedalusMLT", ncFloat, dims);
    NcVar mlt2 = dataFile.addVar("DaedalusMagneticLongitude", ncFloat, dims);
    NcVar mlt3 = dataFile.addVar("DaedalusMagneticLatitude", ncFloat, dims);
    lat.putVar(outlat);
    lat.putAtt("UNITS","deg [-90,90]");
    time.putVar(outtime);
    time.putAtt("UNITS", "Seconds Since 1 Jan 2015 00:00:00.000");
    lon.putVar(outlon);
    lon.putAtt("UNITS","deg [0,180)");
    alt.putVar(outalt);
    alt.putAtt("UNITS","km");
    mlt1.putVar(outacc1);
    mlt1.putAtt("UNITS","magnetic local time");
    mlt2.putVar(outacc2);
    mlt2.putAtt("UNITS","deg");
    mlt3.putVar(outacc3);
    mlt3.putAtt("UNITS","deg");

    delete[] outtime;
    delete[] outlat;
    delete[] outlon;
    delete[] outalt;
    delete[] outacc1;
    delete[] outacc2;
    delete[] outacc3;
    return true;

}



    bool WriteData(std::vector<float> data,std::string varName)
{

    auto out =new  float  [data.size()];
    


    for (int i=0; i<data.size();++i){

        out[i]=data.at(i);
        if (data.at(i)<RUN::fillValue){ out[i]= std::numeric_limits<double>::quiet_NaN();}

    }

    int Size=data.size();
    

    NcFile dataFile(orbit::OutFile, NcFile::write);
    int cnt=dataFile.getDimCount();
    std::cout<<cnt<<std::endl;

    
    vector<NcDim> dim;
    NcDim dims= dataFile.getDim("time");
    
    dim.push_back(dims);


    NcVar outVar = dataFile.addVar(varName,ncFloat,dim);
    outVar.putVar(out);

    if (varName=="BX_si")       {  outVar.putAtt("UNITS","Tesla");}
    if (varName=="BX_si")       {  outVar.putAtt("Long Name","Eastward Component of Magnetic Field");}
    if (varName=="BY_si")       {  outVar.putAtt("UNITS","Tesla");}
    if (varName=="BX_si")       {  outVar.putAtt("Long Name","Northward Component of Magnetic Field");}
    if (varName=="BZ_si")       {  outVar.putAtt("UNITS","Tesla");}
    if (varName=="BX_si")       {  outVar.putAtt("Long Name","Upward Component of Magnetic Field");}
    if (varName=="EEX_si")      {  outVar.putAtt("UNITS","V/m");}
    if (varName=="EEY_si")      {  outVar.putAtt("UNITS","V/m");}
    if (varName=="EEZ_si")      {  outVar.putAtt("UNITS","V/m");}
    if (varName=="SIGMAP_PED")  {  outVar.putAtt("UNITS","S/m");}
    if (varName=="SIGMA_HAL")  {  outVar.putAtt("UNITS","S/m");}
    if (varName=="UN_si")       {  outVar.putAtt("UNITS","m/s");}
    if (varName=="VN_si")       {  outVar.putAtt("UNITS","m/s");}
    if (varName=="Wi_lev_si")   {  outVar.putAtt("UNITS","m/s");}
    if (varName=="Ohmic")       {  outVar.putAtt("UNITS","W/m^3");}
    if (varName=="ZGMID")       {  outVar.putAtt("UNITS","cm");}
    if (varName=="DEN")       {  outVar.putAtt("UNITS","m-3");}
    if (varName=="Convenction_heating")       {  outVar.putAtt("UNITS","W/m^3");}
    if (varName=="Wind_heating")       {  outVar.putAtt("UNITS","W/m^3");}

    delete [] out;

    return true;
}



bool PreaLoadingAvailable(std::array<int,3> currentState){
/*
currentState:
0--> Variable number eg 0 means this is the first variable for interpolation
1--> Pass number eg number of file being read
2--> numPasses total number of files to be read
*/
if (currentState.at(2)<=1){
    return false; //preloading not possible with only 1 file
}

if (currentState.at(0)==0 && currentState.at(1)==0){
    return false; // initialization step 
}

if (currentState.at(1)==currentState.at(2)-1){
    return false; //no next file to be read
}

return true;

}
