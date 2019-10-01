#include "hw2.h"
#include "DataMesh.cpp"
#include "ComputeRHS.cpp"

int main() {
  int N = 1600;
  double Courant = 0.1;
  int Nsteps=N, GhostZone=2;
  double value;
  DataMesh<double> U(GhostZone,{N}), dtU(GhostZone,{N});
  DataMesh<bool> GZ(GhostZone,{N});
  GhostZoneMover GZM({N},GhostZone);
  GZ.SetGZMask(GhostZone,{N});
  /*  for (int i=0; i<204; i++)
    cout << GZ.return_element(i) << ", ";
    cout << endl;*/
  ComputeRHS<double> der({N},GhostZone);
  for (int i = GhostZone; i < N+GhostZone; i++){
    value=sin(2*PI*double(i-GhostZone)/(N));
    U.SetValue(i, value);
  }
  GZM.PeriodicGZ(GZ, U);
  cout <<"{";
  /*  for (int i=0; i<204; i++)
    cout << U.return_element(i) << ", ";
    cout << "}"<<endl;*/
  cout.precision(16);
  for (int j=0; j<=Nsteps; j++) {
    if(j%(Nsteps)==0) {
      cout << "Step=" << j << endl;
      cout << "{";
      for (int i=0; i<N+3;i++)
	cout << U.return_element(i) <<", ";
      cout << U.return_element(N+3) <<"}" << endl;
    }
    der.UpstreamDerivative(U, Courant, dtU,GZ,GZM);
    //der.DownstreamDerivative(U, Courant, dtU,GZ,GZM);
    //der.CenteredDerivative(U, Courant, dtU,GZ,GZM);
    //der.RungeKutta3(U, Courant, dtU, GZ, GZM);
    U+=dtU;
    GZM.PeriodicGZ(GZ, U);
  }
  return 0;
}
