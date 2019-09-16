#include "hw3.h"
#include "DataMesh.cpp"
#include "ComputeRHS.cpp"

int main() {
    int N = 200;
    double Courant = 0.002;
    double value;
    DataMesh<double> U({N}, Courant, 40000);
    for (int i = 0; i < N; i++){
        value=sin(2*PI*double(i)/(N));
        U.SetElement(i, value);
    }
    U.EvolveField();
    return 0;
}
