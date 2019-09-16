#include "hw3.h"


template <typename T>
ComputeRHS<T>::ComputeRHS(vector <int> Size, const double Courant) {
    int step;
    Sizes=Size;
    Npnts=1;
    StencilSteps.push_back(1);
    for (int i =0; i<Size.size(); i++){
        step = 1;
        for (int k=0; k<=i; k++)
            step*=Size[k];
        StencilSteps.push_back(step);
        Npnts*=Size[i];
    }
    for (int i=0; i<Npnts; i++)
        derivative.push_back(0);
    dt=Courant;
}

template <typename T>
ComputeRHS<T>::ComputeRHS(vector <int> Size) {
    int step;
    Npnts=1;
    Sizes=Size;
    StencilSteps.push_back(1);
    for (int i =0; i<Size.size(); i++){
        step = 1;
        for (int k=0; k<=i; k++)
            step*=Size[k];
        StencilSteps.push_back(step);
        Npnts*=Size[i];
    }
    for (int i=0; i<Npnts; i++)
        derivative.push_back(0);
    dt=0.5;
}

template <typename T>
vector<T> ComputeRHS<T>::CenteredDerivative1st(const vector<T>& field) {
    for (int i=0; i<Npnts; i++) {
        derivative[i] = dt*(field[(i + StencilSteps[0]+Sizes[0]) % Sizes[0]] - field[(i - StencilSteps[0]+Sizes[0]) % Sizes[0]]) / 2;
     }
    return derivative;
}

template <typename T>
vector<T> ComputeRHS<T>::DownstreamDerivative1st(const vector<T>& field){
    for (int i=0; i<Npnts; i++)
        derivative[i] = - dt * (field[(i + StencilSteps[0]+Sizes[0]) % Sizes[0]] - field[i % Sizes[0]]);
    return derivative;
}

template <typename T>
vector<T> ComputeRHS<T>::UpstreamDerivative1st(const vector<T>& field){
    for (int i=0; i<Npnts; i++)
        derivative[i] = - dt * (field[i % Sizes[0]] - field[(i - StencilSteps[0]+Sizes[0]) % Sizes[0]]);
    return derivative;
}