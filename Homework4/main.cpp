#include "hw.h"
#include "DataMesh.cpp"
#include "dtU.cpp"

int main() {
  int N = 1601;
  double Courant = 0.1;
  int Nsteps=3*(N-1), GhostZone=2;
  double value;
  vector<int> Size={N};
  DataMesh<double> U(GhostZone,Size), dUt(GhostZone,Size);
  Patch space(Size, {0,1}, GhostZone);
  DataMesh<bool> GZ(GhostZone,Size);
  vector<DataMesh<bool>> FaceGZ(0); 
  GZ.SetGZMask(GhostZone,Size);
  double dx=space.GetCoord(1)[0]-space.GetCoord(0)[0];
  double dt=dx*Courant;
  cout << dt*double(Nsteps) << endl;
  for (int i=0; i<Size.size(); i++){
    DataMesh<bool> FaceGZproj(GhostZone,Size);
    FaceGZproj.SetFaceGZ(GZ,i);
    FaceGZ.push_back(FaceGZproj);
  }
  
  GhostZoneMover GZM(Size,GhostZone);
  dtU<double> der(Size,GhostZone);
  for (int i = GhostZone; i < N+GhostZone; i++){
    value=exp(-(space.GetCoord(i)[0]-0.5)*(space.GetCoord(i)[0]-0.5)/0.04);
    U.SetValue(i, value);
  }
  GZM.GeneratePeriodicGZ(GZ, U);
  GZM.ApplyBCs(U);
  cout.precision(15);
  for (int j=0; j<=Nsteps; j++) {
    if(j%(Nsteps)==0) {
      cout << "Step=" << j << endl;
      cout << "{";
      for (int i=0; i<N+3;i++)
	cout << "{" << space.GetCoord(i)[0] <<", " << U.return_element(i) <<"}, ";
      cout << "{" << space.GetCoord(N+3)[0] <<", " <<U.return_element(N+3) <<"}}" << endl;
    }
    der.RungeKutta3shock(U, Courant, dUt, GZ, GZM, FaceGZ);
    U+=dUt;
    GZM.ApplyBCs(U);
  }
  return 0;
}
