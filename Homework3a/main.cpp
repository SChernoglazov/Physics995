#include "hw2.h"
#include "DataMesh.cpp"
#include "dtU.cpp"

int main() {
  int N = 200;
  double Courant = 0.1;
  int Nsteps=N, GhostZone=2;
  double value;
  DataMesh<double> U(GhostZone,{N}), dUt(GhostZone,{N});
  DataMesh<bool> GZ(GhostZone,{N});
  GhostZoneMover GZM({N},GhostZone);
  GZ.SetGZMask(GhostZone,{N});
  dtU<double> der({N},GhostZone);
  for (int i = GhostZone; i < N+GhostZone; i++){
    value=sin(2*PI*double(i-GhostZone)/(N));
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
    //der.UpstreamDerivative(U, Courant, dUt,GZ,GZM);
    //der.DownstreamDerivative(U, Courant, dUt,GZ,GZM);
    //der.CenteredDerivative(U, Courant, dUt,GZ,GZM);
    der.RungeKutta3(U, Courant, dUt, GZ, GZM);
    U+=dUt;
    GZM.ApplyBCs(U);
  }
  return 0;
}
