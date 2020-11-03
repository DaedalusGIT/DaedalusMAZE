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

#ifndef IO_H
#define IO_H

#include "common.h"

bool ReadGrid(std::vector<float> &glatV,
              std::vector<float> &glonV,
              std::vector<float> &glevV,
              char *sampleFile

);

bool ReadVar(std::vector<float> &varVector,
             std::vector<float> &altVector,
             std::vector<float> &gtimeV,
             char *files[],
             std::string varName,
             std::vector<int> numFile);

bool ReadAsyncVar(std::vector<float> &varVector,
                  std::vector<float> &altVector,
                  std::vector<float> &gtimeV,
                  char *files[],
                  std::string varName,
                  int FileID,
                  bool filesExist);

bool ReadMultiVar(std::vector<float> &varVector,
             std::vector<float> &altVector,
             std::vector<float> &gtimeV,
             char *files[],
             std::string varName,
             std::vector<int> numFile);

bool ReadOrbit(std::vector<float> &lonV,
               std::vector<float> &latV,
               std::vector<float> &altV,
               std::vector<float> &acc1,
               std::vector<float> &acc2,
               std::vector<float> &acc3);

bool WriteOrbit(std::vector<float> dlat,
                std::vector<float> dlon,
                std::vector<float> dalt,
                std::vector<float> dtime,
                std::vector<float> acc1,
                std::vector<float> acc2,
                std::vector<float> acc3);

bool WriteData(std::vector<float> data,
               std::string varName
                );

bool PreaLoadingAvailable(std::array<int,3> currentState);
#endif
