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
#include "common.h"
#include<fstream>

bool logInit(char *files[],int numFiles)
{
    std::ofstream logger;
    logger.open("logfile.txt");
    logger<< "============Interpolator Initialized================" <<std::endl;
    logger<< "MODE:"<<RUN::DEVICE<<std::endl;
    logger<< "METHOD:"<<RUN::METHOD<<std::endl;
    #ifdef PRELOAD_FILES
        logger<< "WARNING: PRELOADING FILES IS ENABLED. USE WITH CAUTION :)";
    #endif
    logger<< "IO Stats:" <<std::endl;
    logger<<"\t"<<"PropagateTime=" <<orbit::PropagateTime<<std::endl;
    logger<<"\t"<< "Orbit Size=" <<orbit::DIM<<std::endl;
    logger<<"\t"<< "Orbit DT=" <<orbit::DT<<"sec"<<std::endl;
    logger<<"\t"<< "Orbit Batch Size=" <<orbit::BATCH<<std::endl;
    logger<<"\t"<< "Model DT=" <<model::DT<<"sec"<<std::endl;
    logger<<"\t"<< "Model Batch Size=" <<model::BATCH<<std::endl;
    logger<<"\t"<< "Grid Of Each File Time x Lev x Lat x Lon=" <<grid::NTIMES<<"x"<<grid::NALT<<"x"<<grid::NLAT<<"x"<<grid::NLON<<std::endl;
    logger << "Orbit File:"<< std::endl;
    logger << "\t" <<orbit::OrbitName<< std::endl;
    logger<< "TIEGCM Files:" <<numFiles-1<<std::endl;
    for (int i=1;i<numFiles;++i){
        logger<<"\t"<< files[i]<<std::endl;
    }
    logger<< "============Interpolator Initialized================" <<std::endl;
    logger<< "Progress Update:" <<std::endl;
    logger.close();

    return true;
}


int index(int i, int j, int k, int z)
{

    int pos = i + j * grid::NTIMES*model::BATCH + k * grid::NTIMES*model::BATCH * grid::NALT + z * grid::NTIMES*model::BATCH * grid::NALT * grid::NLAT;
    return pos;
}








// int index(int i, int j, int k, int z)
// {

//     int pos = i + j * grid::NTIMES + k * grid::NTIMES * grid::NALT + z * grid::NTIMES * grid::NALT * grid::NLAT;
//     return pos;
// }

