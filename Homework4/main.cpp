#include "hw.h"
#include "DataMesh.cpp"
#include "dtU.cpp"

int main() {
  int N = 200;
  double Courant = 0.1;
  int Nsteps=1000, GhostZone=2;
  double value;
  vector<int> Size={N};
  DataMesh<double> U(GhostZone,Size), dUt(GhostZone,Size);

  DataMesh<bool> GZ(GhostZone,Size);
  vector<DataMesh<bool>> FaceGZ(0); 
  GZ.SetGZMask(GhostZone,Size);
  for (int i=0; i<Size.size(); i++){
    DataMesh<bool> FaceGZproj(GhostZone,Size);
    FaceGZproj.SetFaceGZ(GZ,i);
    FaceGZ.push_back(FaceGZproj);
  }
  
  GhostZoneMover GZM(Size,GhostZone);
  dtU<double> der(Size,GhostZone);
  for (int i = GhostZone; i < N+GhostZone; i++){
    //value=sin(2*PI*double(i-GhostZone)/(N));
    value=exp(-(i-N/2)*(i-N/2)/400.0);
    U.SetValue(i, value);
  }
  GZM.GeneratePeriodicGZ(GZ, U);
  GZM.ApplyBCs(U);
  cout.precision(16);
  for (int j=0; j<=Nsteps; j++) {
    if(j%(Nsteps)==0) {
      cout << "Step=" << j << endl;
      cout << "{";
      for (int i=0; i<N+3;i++)
	cout << U.return_element(i) <<", ";
      cout << U.return_element(N+3) <<"}" << endl;
    }
    der.RungeKutta3shock(U, Courant, dUt, GZ, GZM, FaceGZ);
    U+=dUt;
    GZM.ApplyBCs(U);
  }
  return 0;
}
