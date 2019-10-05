#include "hw.h"
#include "ComputeRHS.cpp"

template <typename T>
dtU<T>::dtU (const vector <int> Size,const int GhostZone){  
  Npnts=1;
  int step;
  dim=Size.size();
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
  vector<T> helLR(2*dim,0), hel(dim,0);
  for (int i=0; i<Npnts; i++){
    FluxLR.push_back(helLR);
    ULR.push_back(helLR);
    Flux.push_back(hel);
  }
}

template <typename T>
void dtU<T>::RungeKutta3 (const DataMesh<T>& U, const double dt, DataMesh<T>& dUt, const DataMesh<bool>& GZ, GhostZoneMover& GZM) {

  double c1=0, a11=0, a12=0, a13=0;
  double c2=1.0/2.0, a21=1.0/2.0, a22=0, a23=0;
  double c3=1, a31=-1, a32=2, a33=0;
  double b1=1.0/6.0, b2=2.0/3.0, b3=1.0/6.0;

  vector <T> helper(Npnts), helper1(Npnts), k1(Npnts), k2(Npnts), k3(Npnts);
  for (int i=0; i<Npnts;i++){
    k1[i]=RHS.ThirdDerivative(U, i, GZ);
    helper[i]=U.return_element(i) + dt*a21*k1[i];
  }
  GZM.ApplyBCs(helper);
  for (int i=0; i<Npnts; i++){
    k2[i]=RHS.ThirdDerivative(helper, i, GZ);
    helper1[i]=U.return_element(i)+dt*(a31*k1[i]+a32*k2[i]);
  }
  GZM.ApplyBCs(helper1);
  for (int i=0; i<Npnts; i++) {
    k3[i] = RHS.ThirdDerivative(helper1, i, GZ);
  }
  GZM.ApplyBCs(k3);
  for (int i=0; i<Npnts; i++){
    dUt.SetValue(i, dt*(b1 * k1[i] + b2 * k2[i] + b3 * k3[i]));
  }
}

template <typename T>
void dtU<T>::RungeKutta3shock (const DataMesh<T>& U, const double dt, DataMesh<T>& dUt, const DataMesh<bool>& GZ,GhostZoneMover& GZM, const vector<DataMesh<bool>>& FaceGZ) {

  double c1=0, a11=0, a12=0, a13=0;
  double c2=1.0/2.0, a21=1.0/2.0, a22=0, a23=0;
  double c3=1, a31=-1, a32=2, a33=0;
  double b1=1.0/6.0, b2=2.0/3.0, b3=1.0/6.0;

  vector <T> helper(Npnts), helper1(Npnts), k1(Npnts), k2(Npnts), k3(Npnts);
  FluxLRReconstructor(FaceGZ,U);
  FluxReconstructor(FaceGZ);
  divFlux(GZ,k1);
  for (int i=0; i<Npnts;i++){
    helper[i]=U.return_element(i) + dt*a21*k1[i];
  }
  GZM.ApplyBCs(helper);
  FluxLRReconstructor(FaceGZ,helper);
  FluxReconstructor(FaceGZ);
  divFlux(GZ,k2);
  for (int i=0; i<Npnts; i++){
    helper1[i]=U.return_element(i)+dt*(a31*k1[i]+a32*k2[i]);
  }
  GZM.ApplyBCs(helper1);
  FluxLRReconstructor(FaceGZ,helper1);
  FluxReconstructor(FaceGZ);
  divFlux(GZ,k3);
  GZM.ApplyBCs(k3);
  for (int i=0; i<Npnts; i++){
    dUt.SetValue(i, dt*(b1 * k1[i] + b2 * k2[i] + b3 * k3[i]));
  }
}

//conventionally 2*dim - flux from left to the right in direction dim
template <typename T>
void dtU<T>::FluxLRReconstructor(const vector<DataMesh<bool>>& FaceGZ, const DataMesh<T>& U){  //FaceGZ must be vector of DataMesh<bool> - for each direction
  T uL, uR; 
  for (int i=0; i<dim; i++){
    for (int j=0; j<Npnts; j++){
      if (FaceGZ[i].return_element(j)==0){
	uL=U.return_element(j-StencilSteps[i])+0.5*MinMod(U.return_element(j-2*StencilSteps[i]),U.return_element(j-StencilSteps[i]),U.return_element(j));
	uR=U.return_element(j)-0.5*MinMod(U.return_element(j-StencilSteps[i]),U.return_element(j),U.return_element(j+StencilSteps[i]));
	ULR[j][2*i]=uL;
	ULR[j][2*i+1]=uR;
	//	FluxLR[j][2*i]=RHS.AdvectionFlux(uL);
	//	FluxLR[j][2*i+1]=RHS.AdvectionFlux(uR);
	FluxLR[j][2*i]=RHS.BurgersFlux(uL);
	FluxLR[j][2*i+1]=RHS.BurgersFlux(uR);
      }
    }
  }
}


template <typename T>
void dtU<T>::FluxLRReconstructor(const vector<DataMesh<bool>>& FaceGZ, const vector<T>& U){  //FaceGZ must be vector of DataMesh<bool> - for each direction                 
  T uL, uR;
  for (int i=0; i<dim; i++){
    for (int j=0; j<Npnts; j++){
      if (FaceGZ[i].return_element(j)==0){
        uL=U[j-StencilSteps[i]]+MinMod(U[j-2*StencilSteps[i]],U[j-StencilSteps[i]],U[j])/2;
        uR=U[j]-MinMod(U[j-StencilSteps[i]],U[j],U[j+StencilSteps[i]])/2;
        ULR[j][2*i]=uL;
        ULR[j][2*i+1]=uR;
	//        FluxLR[j][2*i]=RHS.AdvectionFlux(uL);
	//	FluxLR[j][2*i+1]=RHS.AdvectionFlux(uR);
	FluxLR[j][2*i]=RHS.BurgersFlux(uL);
        FluxLR[j][2*i+1]=RHS.BurgersFlux(uR);
      }
    }
  }
}

template <typename T>
int dtU<T>::sign(const T a)const {
  if (a>0)
    return 1;
  else 
    return -1;
}

template <typename T>
T dtU<T>::MinMod(const T a, const T b, const T c) const{
  T result;
  if((c-b)*(b-a)<=0)
    return 0;
  else
    return sign(b-a)*fmin(fabs(b-a),fabs(c-b));
}

template <typename T>
void dtU<T>::FluxReconstructor(const vector<DataMesh<bool>>& FaceGZ){
  T lambda;
  for (int i=0; i<dim; i++){
    for (int j=0; j<Npnts; j++){
      if (FaceGZ[i].return_element(j)==0){
	lambda=fmax(fabs(ULR[j][2*i]),fabs(ULR[j][2*i+1]));
        Flux[j][i]=(FluxLR[j][2*i]+FluxLR[j][2*i+1])/2-lambda*(ULR[j][2*i]-ULR[j][2*i+1])/2; // LLF
      }
    }
  }
}

template <typename T>
void dtU<T>::divFlux(const DataMesh<bool>& GZ,vector<T>&helper){
  T spat;
  for (int i=0; i<Npnts; i++){
    spat=0;
    if (GZ.return_element(i)==0){
      for (int j=0; j<dim; j++){
	spat+=Flux[i+StencilSteps[j]][j]-Flux[i][j];
      }
      helper[i]=spat;
    }
  }
}
