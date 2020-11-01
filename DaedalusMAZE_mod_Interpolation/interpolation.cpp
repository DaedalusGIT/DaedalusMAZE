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
#include "common.h"
#include <omp.h>
#include"cuInterp.h"
#include "coefArray.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include<array>

using namespace std;

bool GaussianLPF(std::vector<float>& intdata ,std::vector<float> dalt,int numPasses){

   std::vector<float> swap ;
   int halfwidth =  1;

  for (int i = 0 ; i< numPasses; i++){
    swap=intdata;
    for (int j =halfwidth+1 ; j< intdata.size()-halfwidth-1; j++ ){
      bool Apply=true;
      //Discard out of range data poits
      for (int c=-2; c<2; c++){
        Apply= Apply &&  intdata.at(j+c)>RUN::fillValue;
      }

      if (!Apply){
        continue;
      }
      float pixel =0.0;
      for (int kk =0; kk<3; ++kk){
        pixel+= intdata.at(j+kk - halfwidth)*blurKernel[kk];
      }
      swap.at(j) = pixel;  
    }
    intdata=swap;
  }
    return true;
}


bool buildDtime(std::vector<float> &dtime){


    float dt= orbit::DT;
    float time=0;
    for (int i=0; i<orbit::DIM; ++i){

        time +=dt;
        dtime.push_back(time);
    }

    return true;
}


bool propagateTime(int &counter,float &time){

 if (orbit::PropagateTime){
    if (time >= model::DT * model::BATCH){
      counter += 1;
      time = 0.0;
      }
  }
 return true;
}


int local(std::vector<float>arr ,float y){

    for (int i=0; i<arr.size()-1; ++i){

        if (y>=arr.at(i)&& y<arr.at(i+1)){
            return i;
        }
    }

    std::cout<<"No Local Neigbors Found..Boundary Condition ERROR"<<std::endl;
    std::cout << "Aborting at Line:" << __LINE__ << std::endl;
    abort();

    return 5;
}




bool dfdt(std::vector<float>f,std::multimap<int, int> dataMap, std::vector<float>& derivative){
  
  std::fill (derivative.begin(),derivative.end(),0.0);

  for (int i=1; i<f.size()-1; ++i){

    // Get LocalPositions
    int p0=dataMap.find(i)->second;
    int pf=dataMap.find(i+1)->second;
    int pb=dataMap.find(i-1)->second;

    bool areInSameOrbit =  (p0+1==pf && p0-1==pb);

    if (!areInSameOrbit){
      continue;
    }
    derivative.at(i) =abs( f.at(i+1) -2*f.at(i) +f.at(i-1));   
    
    bool isSkipped = (f.at(i)<-1.0e20 || f.at(i-1)<-1.0e20 || f.at(i+1)<-1.0e20 );
    
    if (isSkipped){
      derivative.at(i)=0.;
    }
  
  }

  return true;
}


bool sharpEdges(std::vector<float> dfdt,std::vector<int>& points ,float normVal,float tol){

  for (int i=0; i<dfdt.size(); ++i){

    dfdt.at(i)/=normVal;
    // printf("DEr--> %f\n",normVal);
    bool doAssign = (dfdt.at(i)>= tol && i>50 && i<dfdt.size()-50  && i< 10000) ;
    if (doAssign){
      points.push_back(i);
    }

  }

  return true;
}


bool dfdtFix(std::vector<float>& data,std::vector<int> points){

std::vector<float> swap;
int halfWidth = 1;
int window = 10;


for (int i =0; i<points.size(); ++i){
  // if (i %100==0){std::cout<<i<<std::endl;}

  swap=data;

  for (int j = points.at(i)-window; j<points.at(i)+window; ++j){

    float pixel =0.0;
    for (int kk =0; kk<3; ++kk){
      
      pixel+= data.at(j+kk-halfWidth)*blurKernel[kk];

    }
    swap.at(j) = pixel;  
  }
  data=swap;


}

  return true;
}


bool AdaptiveLPF(std::vector<float>& intdata,std::vector<float> dalt,float srt){

  // Multimap Here
  std::vector<float> data;
  std::multimap<int, int> dataMap; 
  int counter=0;

  // Keep only inRange values
  for (int i=0;i<dalt.size();++i){

      bool isInRange = (dalt.at(i) > model::limits::minALT && dalt.at(i) < model::limits::maxALT);

      if (isInRange){
        data.push_back(intdata.at(i));
        dataMap.insert(std::pair<int,int>(counter,i));
        counter++;
      }
  }
  
  std::vector<float> der = data;
  std::vector<int> edges;
  float tol =0.01;
  std::cout<< "Getting in "<<std::endl;
  dfdt(data,dataMap,der);
  float normVal=*max_element(der.begin(), der.end());
  std::cout<< normVal<<std::endl;
  std::cout<< "Edges"<<std::endl;
  sharpEdges(der,edges,normVal,tol);
  std::cout<< edges.size()<<std::endl;
  int numPasses=0;
  bool applyLPF= edges.size()>0;

  while(applyLPF){

    
    numPasses++;
    dfdt(data,dataMap, der);
    sharpEdges(der,edges,normVal,tol);
    applyLPF =  edges.size()>0;
    if (numPasses >10 || !applyLPF) { break;}
    dfdtFix(data, edges);
    cout << edges.size() << endl;
    edges.resize(0);
  }

  

  // Merge based on Multimap Here
  for (int i = 0; i < data.size(); ++i){
      
      int pos = dataMap.find(i)->second;
      intdata.at(pos)=data.at(i);
  }

  return true;
}





bool Derivative(std::vector<float> &f,
                          int time,
                          int x,
                          int y,
                          int z,
                          float &derivative,
                          int mode
                          ){


  int pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8;
  int pos,posP1,posM1;


  if (mode == INTERPOLATION::DERIVATIVES::DOXC){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x+1,y,z);
    posM1 = index(time,x-1,y,z);

    derivative=(f.at(posP1)-f.at(posM1))/2.0;
    return true;
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOXB){
    //Get Indeces
    pos = index(time,x,y,z);
    posM1 = index(time,x-1,y,z);

    derivative=(f.at(pos)-f.at(posM1));
    return true;
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOXF){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x+1,y,z);

    derivative=(f.at(posP1)-f.at(pos));
    return true;
  }


  if (mode == INTERPOLATION::DERIVATIVES::DOYC){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x,y+1,z);
    posM1 = index(time,x,y-1,z);

    derivative=(f.at(posP1)-f.at(posM1))/2.0;
    return true;
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOYB){
    //Get Indeces
    pos = index(time,x,y,z);
    posM1 = index(time,x,y-1,z);

    derivative=(f.at(pos)-f.at(posM1));
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOYF){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x,y+1,z);

    derivative=(f.at(posP1)-f.at(pos));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOZC){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x,y,z+1);
    posM1 = index(time,x,y,z-1);

    derivative=(f.at(posP1)-f.at(posM1))/2.0;
    return true;
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOZB){
    //Get Indeces
    pos = index(time,x,y,z);
    posM1 = index(time,x,y,z-1);

    derivative=(f.at(pos)-f.at(posM1));
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOZF){
    //Get Indeces
    pos = index(time,x,y,z);
    posP1 = index(time,x,y,z+1);

    derivative=(f.at(posP1)-f.at(pos));
    return true;
  }


  if (mode == INTERPOLATION::DERIVATIVES::DOXYC){
    //Get Indeces
    pos1 = index(time,x+1,y+1,z);
    pos2 = index(time,x+1,y-1,z);
    pos3 = index(time,x-1,y+1,z);
    pos4 = index(time,x-1,y-1,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4))/4.0;
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXYF){
    //Get Indeces
    pos1 = index(time,x+1,y+1,z);
    pos2 = index(time,x+1,y,z);
    pos3 = index(time,x,y+1,z);
    pos4 = index(time,x,y,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXYB){
    //Get Indeces
    pos1 = index(time,x,y,z);
    pos2 = index(time,x,y-1,z);
    pos3 = index(time,x-1,y,z);
    pos4 = index(time,x-1,y-1,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXZC){
    //Get Indeces
    pos1 = index(time,x+1,y,z+1);
    pos2 = index(time,x+1,y,z-1);
    pos3 = index(time,x-1,y,z+1);
    pos4 = index(time,x-1,y,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4))/4.0;
    return true;
  }
  if (mode == INTERPOLATION::DERIVATIVES::DOXZF){
    //Get Indeces
    pos1 = index(time,x+1,y,z+1);
    pos2 = index(time,x+1,y,z);
    pos3 = index(time,x,y,z+1);
    pos4 = index(time,x,y,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXZB){
    //Get Indeces
    pos1 = index(time,x,y,z);
    pos2 = index(time,x,y,z-1);
    pos3 = index(time,x-1,y,z);
    pos4 = index(time,x-1,y,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }


  if (mode == INTERPOLATION::DERIVATIVES::DOYZC){
    //Get Indeces
    pos1 = index(time,x,y+1,z+1);
    pos2 = index(time,x,y+1,z-1);
    pos3 = index(time,x,y-1,z+1);
    pos4 = index(time,x,y-1,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4))/4.0;
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOYZF){
    //Get Indeces
    pos1 = index(time,x,y+1,z+1);
    pos2 = index(time,x,y+1,z);
    pos3 = index(time,x,y,z+1);
    pos4 = index(time,x,y,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOYZB){
    //Get Indeces
    pos1 = index(time,x,y,z);
    pos2 = index(time,x,y,z-1);
    pos3 = index(time,x,y-1,z);
    pos4 = index(time,x,y-1,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4));
    return true;
  }


  if (mode == INTERPOLATION::DERIVATIVES::DOXYZC){
    //Get Indeces
    pos1 = index(time,x+1,y+1,z+1);
    pos2 = index(time,x+1,y+1,z-1);
    pos3 = index(time,x+1,y-1,z+1);
    pos4 = index(time,x+1,y-1,z-1);
    pos5 = index(time,x-1,y+1,z+1);
    pos6 = index(time,x-1,y+1,z-1);
    pos7 = index(time,x-1,y-1,z+1);
    pos8 = index(time,x-1,y-1,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4)-f.at(pos5)+f.at(pos6)+f.at(pos7)-f.at(pos8))/8.0;
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXYZF){
    //Get Indeces
    pos1 = index(time,x+1,y+1,z+1);
    pos2 = index(time,x+1,y+1,z);
    pos3 = index(time,x+1,y,z+1);
    pos4 = index(time,x+1,y,z);
    pos5 = index(time,x,y+1,z+1);
    pos6 = index(time,x,y+1,z);
    pos7 = index(time,x,y,z+1);
    pos8 = index(time,x,y,z);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4)-f.at(pos5)+f.at(pos6)+f.at(pos7)-f.at(pos8));
    return true;
  }

  if (mode == INTERPOLATION::DERIVATIVES::DOXYZB){
    //Get Indeces
    pos1 = index(time,x,y,z);
    pos2 = index(time,x,y,z-1);
    pos3 = index(time,x,y-1,z);
    pos4 = index(time,x,y-1,z-1);
    pos5 = index(time,x-1,y,z);
    pos6 = index(time,x-1,y,z-1);
    pos7 = index(time,x-1,y-1,z);
    pos8 = index(time,x-1,y-1,z-1);

    derivative=(f.at(pos1)-f.at(pos2)-f.at(pos3)+f.at(pos4)-f.at(pos5)+f.at(pos6)+f.at(pos7)-f.at(pos8));
    return true;
  }
  exit (EXIT_FAILURE);
  abort();

  return false;
}



bool findSpatialNeighbors(std::vector<float> &glat,
                          std::vector<float> &glon,
                          std::vector<float> &alts,
                          std::vector<float> &zg,
                          float lat,
                          float lon,
                          float alt, int &ltheta,
                          int &lphi, int &lrho,
                          int counter)
{
  //Lon
  if (lon >= glon.back()){
    lphi=glon.size()-1;
  }
  else if (lon <= glon[0]){
    lphi = 0;
  }else{
    lphi = local(glon, lon);
  }

  //Lat
  if (lat >= glat.back()){
    ltheta=glat.size() -1;
  }else if (lat <= glat[0] ){
    ltheta= 0;
  }else{
    ltheta = local(glat, lat);
  }

  // Get Alts from ZGMID
  int dataIndex;
  for (int j =0; j<grid::NALT; ++j){
    dataIndex=index(counter,j,ltheta,lphi);
    alts.at(j) = zg.at(dataIndex) / 1.0e5;
  }

  //Sanity Checks
  if (alts.at(0)>1e20 || alt>= alts.back()){
    return false;
  }

  // Get lrho
  lrho=local(alts,alt);

  return true;
}


bool isBoundaryMAX( int localpos , int dimension){

  bool isBndMax = false;
  switch (dimension){

  case model::dimensions::ALT :
    isBndMax = localpos >= model::limits::maxALT ;
    break;

  case model::dimensions::LON :
    isBndMax = localpos >= model::limits::LONBOUNDARYMAX ;
    break;

  case model::dimensions::LAT :
     isBndMax = localpos >= model::limits::LATBOUNDARYMAX ;
    break;

  default:
    isBndMax = false;
    abort();
    break;
  }

  return isBndMax;
}


bool isBoundaryMIN( int localpos , int dimension){

  bool isBndMin = false;
  switch (dimension){

  case model::dimensions::ALT :
    isBndMin = localpos <=model::limits::minALT  ;
    break;

  case model::dimensions::LON :
    isBndMin = localpos <= model::limits::LONBOUNDARYMIN ;
    break;

  case model::dimensions::LAT :
    isBndMin = localpos <=model::limits::LATBOUNDARYMIN ;
    break;

  default:
    isBndMin = false;
    abort();
    break;
  }

  return isBndMin;
}


bool InterTrilinearOver(int counter,
                        float dx,
                        float dy,
                        float dz,
                        int lrho,
                        int ltheta,
                        int lphi,
                        std::vector<float> &var,
                        std::vector<float> &intData){
    
    std::array<float,8> w;
    float value;
    int dataIndex=0;

    w.at(0)= (1 - dx) * (1 - dy) * (1 - dz);
    w.at(1) = (dx) * (1 - dy) * (1 - dz);
    w.at(2) = (1 - dx) * (dy) * (1 - dz);
    w.at(3) = (dx) * (dy) * (1 - dz);
    w.at(4) = (1 - dx) * (1 - dy) * (dz);
    w.at(5) = (dx) * (1 - dy) * (dz);
    w.at(6) = (1 - dx) * (dy) * (dz);
    w.at(7) = (dx) * (dy) * (dz);

    value  =0.0;
    dataIndex=index(counter,lrho,ltheta,lphi);
    value += w.at(0) * var.at(dataIndex);
    dataIndex=index(counter,lrho+1,ltheta,lphi);
    value += w.at(1) * var.at(dataIndex);
    dataIndex=index(counter,lrho,ltheta+1,lphi);
    value += w.at(2) * var.at(dataIndex);
    dataIndex=index(counter,lrho+1,ltheta+1,lphi);
    value += w.at(3) * var.at(dataIndex);
    dataIndex=index(counter,lrho,ltheta,lphi+1);
    value += w.at(4) * var.at(dataIndex);
    dataIndex=index(counter,lrho+1,ltheta,lphi+1);
    value += w.at(5) * var.at(dataIndex);
    dataIndex=index(counter,lrho,ltheta+1,lphi+1);
    value += w.at(6) * var.at(dataIndex);
    dataIndex=index(counter,lrho+1,ltheta+1,lphi+1);
    value += w.at(7) * var.at(dataIndex);
    intData.push_back(value);
    return true;
}



bool getNeighborDerivativeMode(std::array<int, 8> &modes, int i, int j, int k, int derx, int dery, int derz){


  if (derx ==1 && dery ==0 && derz ==0){ //altitudinal derivatives

    bool notBoundary= !isBoundaryMAX(i+1,model::dimensions::ALT) && !isBoundaryMIN(i-1,model::dimensions::ALT) ;
    bool MaxBoundary= isBoundaryMAX(i+1,model::dimensions::ALT) ;
    bool MinBoundary= isBoundaryMIN(i-1,model::dimensions::ALT) ;

    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOXC);
      return true;
    }

    if (MaxBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOXB;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOXB;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOXB;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOXB;
      return true;
    }

    if (MinBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOXF;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOXF;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOXF;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOXC;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOXF;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOXC;
      return true;
    }
  }

  if (derx ==0 && dery ==1 && derz ==0){ //latitudinal derivatives

    bool notBoundary= !isBoundaryMAX(j+1,model::dimensions::LAT) && !isBoundaryMIN(j-1,model::dimensions::LAT) ;
    bool MaxBoundary= isBoundaryMAX(j+1,model::dimensions::LAT) ;
    bool MinBoundary= isBoundaryMIN(j-1,model::dimensions::LAT) ;

    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOYC);
      return true;
    }

    if (MaxBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOYB;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOYB;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOYB;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOYB;
      return true;
    }

    if (MinBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOYF;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOYF;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOYF;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOYF;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOYC;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOYC;
      return true;
    }
    
  }

  if (derx ==0 && dery ==0 && derz ==1){ //longitudinal derivatives

    bool notBoundary= !isBoundaryMAX(k+1,model::dimensions::LON) && !isBoundaryMIN(k-1,model::dimensions::LON) ;
    bool MaxBoundary= isBoundaryMAX(k+1,model::dimensions::LON) ;
    bool MinBoundary= isBoundaryMIN(k-1,model::dimensions::LON) ;

    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOYC);
      return true;
    }

    if (MaxBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOZB;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOZB;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOZB;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOZB;
      return true;
    }

    if (MinBoundary) {
      modes.at(0)=INTERPOLATION::DERIVATIVES::DOZF;
      modes.at(1)=INTERPOLATION::DERIVATIVES::DOZF;
      modes.at(2)=INTERPOLATION::DERIVATIVES::DOZF;
      modes.at(3)=INTERPOLATION::DERIVATIVES::DOZF;
      modes.at(4)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(5)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(6)=INTERPOLATION::DERIVATIVES::DOZC;
      modes.at(7)=INTERPOLATION::DERIVATIVES::DOZC;
      return true;
    }
  }

  if (derx ==1 && dery ==1 && derz ==0){
    
    bool notBoundary= !isBoundaryMAX(i+2,model::dimensions::ALT) && !isBoundaryMIN(i-2,model::dimensions::ALT)\
                      && !isBoundaryMAX(j+2,model::dimensions::LAT) && !isBoundaryMIN(j-2,model::dimensions::LAT);
    bool MaxBoundary= isBoundaryMAX(i+2,model::dimensions::ALT)||isBoundaryMAX(j+2,model::dimensions::LAT);
    bool MinBoundary= isBoundaryMIN(i-2,model::dimensions::ALT)||isBoundaryMIN(j+2,model::dimensions::LAT);

    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOXYC);
      return true;
    }
    if (MaxBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOXYB;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOXYB; 
      return true;
    }

    if (MinBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOXYF;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOXYF; 
      return true;
    }
    
    
  }
  if (derx ==1 && dery ==0 && derz ==1){
    bool notBoundary= !isBoundaryMAX(i+2,model::dimensions::ALT) && !isBoundaryMIN(i-2,model::dimensions::ALT)\
                      && !isBoundaryMAX(k+2,model::dimensions::LON) && !isBoundaryMIN(k-2,model::dimensions::LON);
    bool MaxBoundary= isBoundaryMAX(i+2,model::dimensions::ALT)||isBoundaryMAX(k+2,model::dimensions::LON);
    bool MinBoundary= isBoundaryMIN(i-2,model::dimensions::ALT)||isBoundaryMIN(k+2,model::dimensions::LON);

    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOXZC);

    return true;
    }

    if (MaxBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOXZB;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOXZB;
    
    return true;
    }

    if (MinBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOXZF;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOXZF;
     return true;
    }
  }

  if (derx ==0 && dery ==1 && derz ==1){
    bool notBoundary= !isBoundaryMAX(j+2,model::dimensions::LAT) && !isBoundaryMIN(j-2,model::dimensions::LAT)\
                      && !isBoundaryMAX(k+2,model::dimensions::LON) && !isBoundaryMIN(k-2,model::dimensions::LON);
    bool MaxBoundary= isBoundaryMAX(j+2,model::dimensions::LAT)||isBoundaryMAX(k+2,model::dimensions::LON);
    bool MinBoundary= isBoundaryMIN(j-2,model::dimensions::LAT)||isBoundaryMIN(k+2,model::dimensions::LON);
    
    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOYZC);
      return true;   
    }

    if (MaxBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOYZB;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOYZB;
      return true;   
    }

    if (MinBoundary) {
      modes.at(0) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(1) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(2) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(3) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(4) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(5) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(6) = INTERPOLATION::DERIVATIVES::DOYZF;
      modes.at(7) = INTERPOLATION::DERIVATIVES::DOYZF;
      return true;   
    }
  }
  if (derx ==1 && dery ==1 && derz ==1){
    bool notBoundary= !isBoundaryMAX(i+3,model::dimensions::ALT) && !isBoundaryMIN(i-3,model::dimensions::ALT)\
                      &&!isBoundaryMAX(j+3,model::dimensions::LAT) && !isBoundaryMIN(j-3,model::dimensions::LAT)\
                      && !isBoundaryMAX(k+3,model::dimensions::LON) && !isBoundaryMIN(k-3,model::dimensions::LON);
    if (notBoundary) {
      modes.fill(INTERPOLATION::DERIVATIVES::DOXYZC);
      return true;
    }
  }

  return false;
}



bool offsetNeighborsExist(std::array<int,3> position,int offset){


  //Altitudinal Check
  bool altBndMIN=  isBoundaryMIN(position.at(0)-offset,model::dimensions::ALT);
  bool altBndMAX=  isBoundaryMAX(position.at(0)+offset,model::dimensions::ALT);


  //Latitudinal Check
  bool thetaBndMIN=  isBoundaryMIN(position.at(1)-offset,model::dimensions::LAT);
  bool thetaBndMAX=  isBoundaryMAX(position.at(1)+offset,model::dimensions::LAT);

  //Longitudinal Check
  bool phiBndMIN=  isBoundaryMIN(position.at(2)-offset,model::dimensions::LON);
  bool phiBndMAX=  isBoundaryMAX(position.at(2)+offset,model::dimensions::LON);

  bool neighborsExist = !(altBndMAX || altBndMIN || thetaBndMAX || thetaBndMIN || phiBndMAX || phiBndMIN);  

  return neighborsExist;
}






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

)
{
    std::vector<float> alts;
    std::array<float,64> coefs,weights;
    std::array<float,8> fval,dfdx,dfdy,dfdz,dfdxdy,dfdxdz,dfdydz,dfdxdydz;
    alts.resize(grid::NALT);
    float time = 0.0,value, dx, dy, dz, deltarho, deltaphi = 2.5, deltatheta = 2.5;
    int counter=0;
    //Open Error Log
    std::ofstream errorFile;
    errorFile.open("errorFile.txt",std::ios_base::app);


    for (int i=0;i<dlat.size();++i){

      bool isInRange = (dalt.at(i) >110.  && dalt.at(i) < 450.);

      if (isInRange){

            int ltheta,lphi,lrho;

            /*Find Spatial Neihbors*/
            if (!findSpatialNeighbors(glat,glon,alts,zg,dlat.at(i),dlon.at(i),dalt.at(i),ltheta,lphi,lrho,counter)){
              errorFile<<"NEIG::\tNeighbors not found for Alt,Lat,Lon,lrho,ltheta,lrho:\t"<<dalt[i]<<"\t"<<dlat[i]<<"\t"<<dlon[i]<<"\t"<<"\t"<<lrho<<"\t"<<ltheta<<"\t"<<lphi<<"\t"<<std::endl;
            }
            
            /*Scheme Availability*/
            std::array<int,3> positions={lrho,ltheta,lphi}; 
            bool thirdNeighborsExist= offsetNeighborsExist(positions,3);
            bool firstNeighborsExist= offsetNeighborsExist(positions,1);
              
	    //std::cout<<alts.at(0)<<std::endl;
	    //std::cout<< ltheta<<"\t"<<lphi<<"\t"<<lrho<<std::endl;
	    //std::cout<<dalt.at(i)<<std::endl;
          if (alts.at(0) > 1.e20){

              errorFile<<"ERROR:: Possible uninitialized Values..." <<std::endl;
             goto SKIP;
      }		

            /*Relative Distances*/
            deltarho=alts.at(lrho+1)-alts.at(lrho);
            dx =((dalt.at(i) - alts.at(lrho)) / deltarho);
            dy =((dlat.at(i) - glat.at(ltheta)) / deltatheta);
            dz =((dlon.at(i) - glon.at(lphi)) / deltaphi);


            if (!thirdNeighborsExist && !firstNeighborsExist){
              errorFile<<"DER::\t 1st and third Neighbors not found  for Alt,Lat,Lon,lrho,ltheta,lrho:\t"<<dalt[i]<<"\t"<<dlat[i]<<"\t"<<dlon[i]<<"\t"<<"\t"<<lrho<<"\t"<<ltheta<<"\t"<<lphi<<"\t"<<std::endl;
              goto SKIP;

            }else if(thirdNeighborsExist==false && firstNeighborsExist==true){
              errorFile<<"TRI::\t Trilinear used for Alt,Lat,Lon,lrho,ltheta,lrho:\t"<<dalt[i]<<"\t"<<dlat[i]<<"\t"<<dlon[i]<<"\t"<<"\t"<<lrho<<"\t"<<ltheta<<"\t"<<lphi<<"\t"<<std::endl;
              InterTrilinearOver(counter,abs(dx),abs(dy),abs(dz),lrho,ltheta,lphi,var,intdata);
              time += orbit::DT;
              propagateTime(counter,time);
              continue;

            } 


            fval.at(0) = var.at(index(counter, lrho, ltheta, lphi));
            fval.at(1)=var.at(index(counter,lrho+1,ltheta,lphi));
            fval.at(2)=var.at(index(counter,lrho,ltheta+1,lphi));
            fval.at(3)=var.at(index(counter,lrho+1,ltheta+1,lphi));
            fval.at(4)=var.at(index(counter,lrho,ltheta,lphi+1));
            fval.at(5)=var.at(index(counter,lrho+1,ltheta,lphi+1));
            fval.at(6)=var.at(index(counter,lrho,ltheta+1,lphi+1));
            fval.at(7)=var.at(index(counter,lrho+1,ltheta+1,lphi+1));

            std::array<int,8> modes;

            modes.fill(INTERPOLATION::DERIVATIVES::DOXC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdx.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdx.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi, dfdx.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdx.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdx.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdx.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdx.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdx.at(7),modes.at(7));

            modes.fill(INTERPOLATION::DERIVATIVES::DOYC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdy.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdy.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdy.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdy.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdy.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdy.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdy.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdy.at(7),modes.at(7));

            modes.fill(INTERPOLATION::DERIVATIVES::DOZC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdz.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdz.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdz.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdz.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdz.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdz.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdz.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdz.at(7),modes.at(7));

            modes.fill(INTERPOLATION::DERIVATIVES::DOXYC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdxdy.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdxdy.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdxdy.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdxdy.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdxdy.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdxdy.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdxdy.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdxdy.at(7),modes.at(7));

            modes.fill(INTERPOLATION::DERIVATIVES::DOXZC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdxdz.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdxdz.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdxdz.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdxdz.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdxdz.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdxdz.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdxdz.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdxdz.at(7),modes.at(7));

            modes.fill(INTERPOLATION::DERIVATIVES::DOYZC);
            Derivative(var,counter,lrho,ltheta,lphi,dfdydz.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdydz.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdydz.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdydz.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdydz.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdydz.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdydz.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdydz.at(7),modes.at(7));

            Derivative(var,counter,lrho,ltheta,lphi,dfdxdydz.at(0),modes.at(0));
            Derivative(var,counter,lrho+1,ltheta,lphi,dfdxdydz.at(1),modes.at(1));
            Derivative(var,counter,lrho,ltheta+1,lphi,dfdxdydz.at(2),modes.at(2));
            Derivative(var,counter,lrho+1,ltheta+1,lphi,dfdxdydz.at(3),modes.at(3));
            Derivative(var,counter,lrho,ltheta,lphi+1,dfdxdydz.at(4),modes.at(4));
            Derivative(var,counter,lrho+1,ltheta,lphi+1,dfdxdydz.at(5),modes.at(5));
            Derivative(var,counter,lrho,ltheta+1,lphi+1,dfdxdydz.at(6),modes.at(6));
            Derivative(var,counter,lrho+1,ltheta+1,lphi+1,dfdxdydz.at(7),modes.at(7));

            for (int ii=0; ii<8; ++ii){
                coefs.at(ii) = fval.at(ii);
                coefs.at(ii+8) = dfdx.at(ii);
                coefs.at(ii+16) = dfdy.at(ii);
                coefs.at(ii+24) = dfdz.at(ii);
                coefs.at(ii+32) = dfdxdy.at(ii);
                coefs.at(ii+40) = dfdxdz.at(ii);
                coefs.at(ii+48) = dfdydz.at(ii);
                coefs.at(ii+56) = dfdxdydz.at(ii);

            }


            for (int ii=0; ii<64; ++ii){
                weights.at(ii)=0.0;
                for (int jj=0; jj<64; ++jj){

                    weights.at(ii)+=A[ii][jj]*coefs.at(jj);

                }
             }


            value=0;
            for (int ii=0; ii<4; ++ii){
                for (int jj=0; jj<4; ++jj){
                    for (int kk=0; kk<4; ++kk){

                    value+= weights.at(ii+4*jj+16*kk)*pow(dx,ii)*pow(dy,jj)*pow(dz,kk);

                    }
                }
            }
            
            
            intdata.push_back(value);


            time += orbit::DT;
            propagateTime(counter,time);
      }else{

SKIP:intdata.push_back(RUN::fillValue*1.0e5);
            time += orbit::DT;
            propagateTime(counter,time);


        }
    }


    GaussianLPF(intdata,dalt,7);
    return true;
}




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


)
{
    std::array<float,8> w;
    std::vector<float> alts;
    alts.resize(grid::NALT);
    float time = 0.0,value, dx, dy, dz, deltarho, deltaphi = 2.5, deltatheta = 2.5;
    int dataIndex,counter=0;


    for (int i=0;i<dlat.size();++i){


        if (dalt.at(i) > 110. && dalt.at(i) < 480.){

            int ltheta,lphi,lrho;

            if (dlon[i] >= glon[grid::NLON-1]){
                lphi=142;

            }
            else if (dlon[i] <= glon[0]){
                lphi = 0;
            }else{
                lphi = local(glon, dlon.at(i));
            }


            if (dlat[i] >= glat[grid::NLAT-1]){
                ltheta=70;

            }else if (dlat[i] <= glat[0]){
                ltheta=0;
            }else{
                ltheta = local(glat, dlat.at(i));
            }



            for (int j =0; j<grid::NALT; ++j){

                dataIndex=index(counter,j,ltheta,lphi);
                alts.at(j) = zg.at(dataIndex) / 1.0e5;
                if (alts.at(0)>1e20){
                    goto SKIP;
                }
            }

            if (dalt.at(i)>= alts.at(grid::NALT-1)){
                goto SKIP;
            }


            lrho=local(alts,dalt.at(i));
            deltarho=alts.at(lrho+1)-alts.at(lrho);


            dx = abs((dalt.at(i) - alts.at(lrho)) / deltarho);
            dy = abs((dlat.at(i) - glat.at(ltheta)) / deltatheta);
            dz = abs((dlon.at(i) - glon.at(lphi)) / deltaphi);


            deltarho=alts.at(lrho+1)-alts.at(lrho);
            dx =((dalt.at(i) - alts.at(lrho)) / deltarho);
            dy =((dlat.at(i) - glat.at(ltheta)) / deltatheta);
            dz =((dlon.at(i) - glon.at(lphi)) / deltaphi);

            w.at(0)= (1 - dx) * (1 - dy) * (1 - dz);
            w.at(1) = (dx) * (1 - dy) * (1 - dz);
            w.at(2) = (1 - dx) * (dy) * (1 - dz);
            w.at(3) = (dx) * (dy) * (1 - dz);
            w.at(4) = (1 - dx) * (1 - dy) * (dz);
            w.at(5) = (dx) * (1 - dy) * (dz);
            w.at(6) = (1 - dx) * (dy) * (dz);
            w.at(7) = (dx) * (dy) * (dz);

            value  =0.0;
            dataIndex=index(counter,lrho,ltheta,lphi);
            value += w.at(0) * var.at(dataIndex);
            dataIndex=index(counter,lrho+1,ltheta,lphi);
            value += w.at(1) * var.at(dataIndex);
            dataIndex=index(counter,lrho,ltheta+1,lphi);
            value += w.at(2) * var.at(dataIndex);
            dataIndex=index(counter,lrho+1,ltheta+1,lphi);
            value += w.at(3) * var.at(dataIndex);
            dataIndex=index(counter,lrho,ltheta,lphi+1);
            value += w.at(4) * var.at(dataIndex);
            dataIndex=index(counter,lrho+1,ltheta,lphi+1);
            value += w.at(5) * var.at(dataIndex);
            dataIndex=index(counter,lrho,ltheta+1,lphi+1);
            value += w.at(6) * var.at(dataIndex);
            dataIndex=index(counter,lrho+1,ltheta+1,lphi+1);
            value += w.at(7) * var.at(dataIndex);
            intdata.push_back(value);
            time += orbit::DT;

            if (orbit::PropagateTime){
                if (time >= model::DT * model::BATCH)
                {
                    counter += 1;
                    time = 0.0;
                }
            }
        }else{
            SKIP:intdata.push_back(RUN::fillValue*1.0e5);
            time += orbit::DT;
            if (orbit::PropagateTime)
            {
                if (time >= model::DT*model::BATCH)
                {
                    counter += 1;
                    time = 0.0;
                }
            }
        }
        
    }
    
    return true;
}

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
                      std::vector<float> &intdata)

{


    auto gridLat = new float [grid::NLAT];
    auto gridLon = new float [grid::NLON];
    auto gridLev = new float [grid::NALT];
    auto gridTime = new float [grid::NTIMES*model::BATCH];

    auto daedLat  = new float [orbit::BATCH];
    auto daedVar  = new float [orbit::BATCH];
    auto daedLon  = new float [orbit::BATCH];
    auto daedAlt  = new float [orbit::BATCH];
    auto gridAlts = new float [grid::NTIMES*model::BATCH*grid::NLON*grid::NLAT*grid::NALT];
    auto gridVar  = new float [grid::NTIMES*model::BATCH*grid::NLON*grid::NLAT*grid::NALT];


    for (int i=0; i<orbit::BATCH; ++i){

        daedLon[i]=dlon.at(i);
        daedLat[i]=dlat.at(i);
        daedAlt[i]=dalt.at(i);

    }


    for (int i =0; i<grid::NTIMES*model::BATCH*grid::NLON*grid::NLAT*grid::NALT; ++i){


        if (i<grid::NTIMES*model::BATCH){gridTime[i]=gtime.at(i);}
        if (i<grid::NLAT){gridLat[i]=glat.at(i);}
        if (i<grid::NLON){gridLon[i]=glon.at(i);}
        if (i<grid::NALT){gridLev[i]=glev.at(i);}

        gridAlts[i]=zg.at(i);
        gridVar[i]=var.at(i);

    }

    DeviceInterpolateBatch(gridLat,gridLon,gridLev,gridTime,daedLat,daedLon,daedAlt,gridAlts,gridVar,daedVar);

    for (int i=0; i<orbit::BATCH; ++i){

        intdata.push_back(daedVar[i]);

    }


    delete []  gridLat ;
    delete []  gridLon ;
    delete []  gridLev ;
    delete []  gridTime;
    delete []  daedLat ;
    delete []  daedVar ;
    delete []  daedLon ;
    delete []  daedAlt ;
    delete []  gridAlts;
    delete []  gridVar ;



    return true;
}

#endif
bool Boxcar( float *&data,int Size){


    int filterPasses=2;
    auto swap=new float[Size];




    for (int j=0;j<filterPasses;++j){

        for (int i=1;i<Size-1;++i){

            swap[i]=0.33*(data[i-1]+data[i]+data[i+1]);


        }
        for (int i = 1; i < Size - 1; ++i){
            data[i]=swap[i];
        }
    }


    return true;
  }

