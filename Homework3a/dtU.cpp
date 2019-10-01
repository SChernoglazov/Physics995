#include "hw2.h"
#include "ComputeRHS.cpp"

template <typename T>
dtU<T>::dtU (vector <int> Size, int GhostZone){  
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
  RHS.FillComputeRHS(Size, GhostZone);
}

template <typename T>
void dtU<T>::RungeKutta3 (DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM) {

  double c1=0, a11=0, a12=0, a13=0;
  double c2=1.0/2.0, a21=1.0/2.0, a22=0, a23=0;
  double c3=1, a31=-1, a32=2, a33=0;
  double b1=1.0/6.0, b2=2.0/3.0, b3=1.0/6.0;

  vector <T> helper(Npnts), helper1(Npnts), k1(Npnts), k2(Npnts), k3(Npnts);
  for (int i=0; i<Npnts;i++){
    k1[i]=RHS.ThirdDerivative(U, i, GZ);
    helper[i]=U.return_element(i) + dt*a21*k1[i];
  }
  GZM.PeriodicGZ(GZ, helper);
  for (int i=0; i<Npnts; i++){
    k2[i]=RHS.ThirdDerivative(helper, i, GZ);
    helper1[i]=U.return_element(i)+dt*(a31*k1[i]+a32*k2[i]);
  }
  GZM.PeriodicGZ(GZ,helper1);
  for (int i=0; i<Npnts; i++) {
    k3[i] = RHS.ThirdDerivative(helper1, i, GZ);
  }
  GZM.PeriodicGZ(GZ,k3);
  for (int i=0; i<Npnts; i++){
    dUt.SetValue(i, dt*(b1 * k1[i] + b2 * k2[i] + b3 * k3[i]));
  }
}

template <typename T>
void dtU<T>::UpstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM){
  T current_point, previous_point, result;
  for (int i=0; i<Npnts; i++){
    if(GZ.return_element(i)==0){
      current_point=U.return_element(i);
      previous_point=U.return_element(i-StencilSteps[0]);
      result=-dt*(current_point-previous_point);
      dUt.SetValue(i,result);
    }
  }
  GZM.PeriodicGZ(GZ, dUt);
}

template <typename T>
void dtU<T>::DownstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt,DataMesh<bool>& GZ, GhostZoneMover& GZM){
  T current_point, next_point, result;
  for (int i=0; i<Npnts; i++){
    if(GZ.return_element(i)==0){
      current_point=U.return_element(i);
      next_point=U.return_element(i+StencilSteps[0]);
      result=-dt*(next_point-current_point);
      dUt.SetValue(i,result);
    }
  }
  GZM.PeriodicGZ(GZ, dUt);
}

template <typename T>
void dtU<T>::CenteredDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt,DataMesh<bool>& GZ, GhostZoneMover& GZM){
  T previous_point, next_point, result;
  for (int i=0; i<Npnts; i++){
    if(GZ.return_element(i)==0){
      previous_point=U.return_element(i-StencilSteps[0]);
      next_point=U.return_element(i+StencilSteps[0]);
      result=-dt*(next_point-previous_point)/2;
      dUt.SetValue(i,result);
    }
  }
  GZM.PeriodicGZ(GZ, dUt);
}
