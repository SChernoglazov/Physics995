#include "hw3.h"
#include "DataMesh.cpp"
#include "ComputeRHS.cpp"

int main() {
    int N = 800;
    double Courant = 0.1*1/N;
    int Nsteps=10*N;
    //int Nsteps=1;
    double value;
//the wave crosses the domain 5 times exactly
    DataMesh<double> U({N}), dtU({N});
    ComputeRHS<double> der({N});
    for (int i = 0; i < N; i++){
        value=sin(2*PI*double(i)/(N));
        U.SetElement(i, value);
    }
    for (int j=0; j<=Nsteps; j++) {
      if(j%(Nsteps)==0) {
            cout << "Step=" << j << endl;
            cout << "{";
            for (int i=0; i<N-1;i++)
                cout <<U.GetElement(i) <<", ";
            cout << U.GetElement(N-1) <<"}" << endl;
        }
      //der.UpstreamDerivative(U, Courant, dtU);
        der.DownstreamDerivative(U, Courant, dtU);
        //der.CenteredDerivative(U, Courant, dtU);
        //der.RungeKutta3(U, Courant, dtU);
        U+=dtU;
    }
    return 0;
}
