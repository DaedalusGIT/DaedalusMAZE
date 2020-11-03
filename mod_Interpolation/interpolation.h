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

#ifndef INTERPOLATION_TOOLS_H
#define INTERPOLATION_TOOLS_H

#include "common.h"
#include <map>
bool InterpolateBatch(std::vector<float> &glat,
                      std::vector<float> &glon,
                      std::vector<float> &glev,
                      std::vector<float> &gtime,
                      std::vector<float> &dlat,
                      std::vector<float> &dlon,
                      std::vector<float> &dalt,
                      std::vector<float> &zg,
                      std::vector<float> &var,
                      std::vector<float> &intdata


);

#ifdef GPU
bool PrepeareForGPU(std::vector<float> &glat,
                      std::vector<float> &glon,
                      std::vector<float> &glev,
                      std::vector<float> &gtime,
                      std::vector<float> &dlat,
                      std::vector<float> &dlon,
                      std::vector<float> &dalt,
                      std::vector<float> &zg,
                      std::vector<float> &var,
                      std::vector<float> &intdata);
#endif


bool isBoundaryMAX( int localpos , int dimension);
bool isBoundaryMIN( int localpos , int dimension);
bool dfdt(std::vector<float> f, std::multimap<int, int> dataMap, std::vector<float> &derivative);
bool sharpEdges(std::vector<float> dfdt,std::vector<int>& points ,float normVal,float tol);
bool AdaptiveLPF(std::vector<float> &data, std::vector<float> dalt,float srt);
bool dfdtFix(std::vector<float>& data,std::vector<float> points);

bool getNeighborDerivativeMode(std::array<int,8>& modes, int i, int j, int k, int derx, int dery, int derz);


bool findSpatialNeighbors(std::vector<float> &glat,
                          std::vector<float> &glon,
                          std::vector<float> &alts,
                          std::vector<float> &zg,
                          float lat,
                          float lon,
                          float alt, int &ltheta,
                          int &lphi, int &lrho,
                          int counter);

bool buildDtime(std::vector<float> &dtime);

int local(std::vector<float> arr, float y);

// bool Boxcar(float data[],int Size);
bool Boxcar(float *&data,int Size);

bool InterTrilinearOver(int counter,
                        float dx,
                        float dy,
                        float dz,
                        int lrho,
                        int ltheta,
                        int lphi,
                        std::vector<float> &var,
                        std::vector<float> &intData);



bool Derivative(std::vector<float> &f,
                          int time,
                          int x,
                          int y,
                          int z,
                          float& derivative,
                          int mode
                          );

bool TriCubicSpline(std::vector<float> &glat,
                      std::vector<float> &glon,
                      std::vector<float> &glev,
                      std::vector<float> &gtime,
                      std::vector<float> &dlat,
                      std::vector<float> &dlon,
                      std::vector<float> &dalt,
                      std::vector<float> &zg,
                      std::vector<float> &var,
                      std::vector<float> &intdata,
                      std::string varname


);


bool offsetNeighborsExist(std::array<int,3> position,int offset);


bool propagateTime(int &counter,float &time);

bool GaussianLPF(std::vector<float>& intdata ,std::vector<float> dalt,int numPasses);

#endif
