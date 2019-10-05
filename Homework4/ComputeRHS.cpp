#include "hw.h"

template <typename T>
ComputeRHS<T>::ComputeRHS(){}

template <typename T>
void ComputeRHS<T>::FillComputeRHS(const vector <int> Size, const int GhostZone){
  Npnts=1;
  int step;
  Sizes=Size;
  GhostZones = GhostZone;
  for (int i=0; i<Size.size();i++)
    Npnts *= (Size[i]+2*GhostZone);
  StencilSteps.push_back(1);
  for (int i =0; i<Size.size(); i++){
    step = 1;
    for (int k=0; k<=i; k++)
      step*=(Size[k]+2*GhostZone);
    StencilSteps.push_back(step);
  }
}

template <typename T>
T ComputeRHS<T>::ThirdDerivative(const vector <T>& helper, const int i, const  DataMesh<bool>& GZ){
  T pointm2, pointm1,pointp1,pointp2, result;
  result=0;
  if (GZ.return_element(i)==0){ //if the point is not masked 
    pointm2=helper[i-2*StencilSteps[0]];
    pointm1=helper[i-StencilSteps[0]];
    pointp2=helper[i+2*StencilSteps[0]];
    pointp1=helper[i+StencilSteps[0]];
    result=-(pointm2-8.0*pointm1+8.0*pointp1-pointp2)/12;
  }
  return result;
}

template <typename T>
T ComputeRHS<T>::ThirdDerivative(const DataMesh<T>& U, const int i, const DataMesh<bool>& GZ) {
  T pointm2, pointm1,pointp1,pointp2, result;
  result=0;
  if(GZ.return_element(i)==0){
    pointm2=U.return_element((i-2*StencilSteps[0]+Npnts)%Npnts);
    pointm1=U.return_element((i-StencilSteps[0]+Npnts)%Npnts);
    pointp1=U.return_element((i+StencilSteps[0]+Npnts)%Npnts);
    pointp2=U.return_element((i+2*StencilSteps[0]+Npnts)%Npnts);
    result=-(pointm2-8.0*pointm1+8.0*pointp1-pointp2)/12;
  }
  return result;
}
template <typename T>
T ComputeRHS<T>::BurgersFlux(const T u)const{
  return -u*u/2;
}

template <typename T>
T ComputeRHS<T>::AdvectionFlux(const T u)const{
  return -u;
}
