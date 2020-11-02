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

#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <iostream>
#include<vector>

bool logInit(char *files[],int numFiles);

int index(int i, int j, int k, int z);

namespace RUN{
    const std::string DEVICE="CPU";
     //const std::string METHOD="TRILINEAR";
    const std::string METHOD="TRICUBIC";
    const bool MANUAL=false;
    // const bool MANUAL=false;
    const float fillValue = -1.0e20;


}

namespace orbit{

    const bool PropagateTime=true;
    const int DIM =;//total orbit datapoints
    const float DT=10.0;      //orbit dtd
    const int BATCH =;    //orbit batch processing
    const int SR=0.1;
    const std::string OrbitName =;
    const std::string OutFile = "IntDataOut.nc";
}

namespace grid{
    const int NTIMES =;
    const int NLAT=72;
    const int NLON=144;
    const int NALT=57;

}


namespace model {

    const int DT=7200;
    const int BATCH=1;

    enum limits{
        maxALT=56,
        minALT=1,
        LONBOUNDARYMAX = 143,
        LONBOUNDARYMIN = 0,
        LATBOUNDARYMAX = 71,
        LATBOUNDARYMIN = 0,
        N_LIMITS
    };
    enum dimensions{
        TIME,
        ALT,
        LAT,
        LON,
        N_DIMENSIONS
    };




}


namespace INTERPOLATION{

  enum DERIVATIVES{
  DOXC,
  DOYC,
  DOZC,
  DOXYC,
  DOXZC,
  DOYZC,
  DOXYZC,
  DOXB,
  DOYB,
  DOZB,
  DOXYB,
  DOXZB,
  DOYZB,
  DOXYZB,
  DOXF,
  DOYF,
  DOZF,
  DOXYF,
  DOXZF,
  DOYZF,
  DOXYZF,
  DOTIME,
  N_DERIVATIVES

  };

}

#endif
