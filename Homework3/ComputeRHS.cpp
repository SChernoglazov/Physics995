#include "hw3.h"

template <typename T>
ComputeRHS<T>::ComputeRHS(vector <int> Size){
    Npnts=1;
    int step;
    for (int i=0; i<Size.size();i++)
        Npnts *= Size[i];
    for (int i=0;i<Npnts;i++)
        derivative.push_back(0);
    StencilSteps.push_back(1);
    for (int i =0; i<Size.size(); i++){
        step = 1;
        for (int k=0; k<=i; k++)
            step*=Size[k];
        StencilSteps.push_back(step);
    }
}

template <typename T>
void ComputeRHS<T>::UpstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dtU){
    T current_point, previous_point, result;
    for (int i=0; i<Npnts; i++){
        current_point=U.GetElement(i);
        previous_point=U.GetElement((i-StencilSteps[0]+Npnts)%Npnts);
        result=-dt*(current_point-previous_point);
        dtU.SetElement(i,result);
    }
}

template <typename T>
void ComputeRHS<T>::DownstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dtU){
    T current_point, next_point, result;
    for (int i=0; i<Npnts; i++){
        current_point=U.GetElement(i);
        next_point=U.GetElement((i+StencilSteps[0]+Npnts)%Npnts);
        result=-dt*(next_point-current_point);
        dtU.SetElement(i,result);
    }
}

template <typename T>
void ComputeRHS<T>::CenteredDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dtU){
    T previous_point, next_point, result;
    for (int i=0; i<Npnts; i++){
        previous_point=U.GetElement((i-StencilSteps[0]+Npnts)%Npnts);
        next_point=U.GetElement((i+StencilSteps[0]+Npnts)%Npnts);
        result=-dt*(next_point-previous_point)/2;
        dtU.SetElement(i,result);
    }
}

template <typename T>
void ComputeRHS<T>::RungeKutta3 (DataMesh<T>& U, const double dt, DataMesh<T>& dtU) {
  /*    double c1=0, a11=0, a12=0, a13=0;
    double c2=1.0/3.0, a21=1.0/3.0, a22=0, a23=0;
    double c3=2.0/3.0, a31=0, a32=2.0/3.0, a33=0;
    double b1=1.0/4.0, b2=0, b3=3.0/4.0;*/

    double c1=0, a11=0, a12=0, a13=0;
    double c2=1.0/2.0, a21=1.0/2.0, a22=0, a23=0;
    double c3=1, a31=-1, a32=2, a33=0;
    double b1=1.0/6.0, b2=2.0/3.0, b3=1.0/6.0;

    vector <T> helper(Npnts), helper1(Npnts), k1(Npnts), k2(Npnts), k3(Npnts);
    for (int i=0; i<Npnts;i++){
        k1[i]=ThirdDerivative(U, i);
        helper[i]=U.GetElement(i) + dt*a21*k1[i];
    }
    for (int i=0; i<Npnts; i++){
        k2[i]=ThirdDerivative(helper, i);
        helper1[i]=U.GetElement(i)+dt*(a31*k1[i]+a32*k2[i]);
    }
    for (int i=0; i<Npnts; i++) {
        k3[i] = ThirdDerivative(helper1, i);
        dtU.SetElement(i, dt*(b1 * k1[i] + b2 * k2[i] + b3 * k3[i]));
    }
}

template <typename T>
T ComputeRHS<T>::ThirdDerivative(vector <T>& helper, const int i){
    T pointm2, pointm1,pointp1,pointp2, result;
    pointm2=helper[(i-2*StencilSteps[0]+Npnts)%Npnts];
    pointm1=helper[(i-StencilSteps[0]+Npnts)%Npnts];
    pointp2=helper[(i+2*StencilSteps[0]+Npnts)%Npnts];
    pointp1=helper[(i+StencilSteps[0]+Npnts)%Npnts];
    result=-(pointm2-8.0*pointm1+8.0*pointp1-pointp2)/12;
    return result;
}

template <typename T>
T ComputeRHS<T>::ThirdDerivative(DataMesh<T>& U, const int i) {
    T pointm2, pointm1,pointp1,pointp2, result;
    pointm2=U.GetElement((i-2*StencilSteps[0]+Npnts)%Npnts);
    pointm1=U.GetElement((i-StencilSteps[0]+Npnts)%Npnts);
    pointp1=U.GetElement((i+StencilSteps[0]+Npnts)%Npnts);
    pointp2=U.GetElement((i+2*StencilSteps[0]+Npnts)%Npnts);
    result=-(pointm2-8.0*pointm1+8.0*pointp1-pointp2)/12;
    return result;
}
