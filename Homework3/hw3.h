#ifndef HW3_HW3_H
#define HW3_HW3_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <math.h>

#define PI 3.141592654
using namespace std;

class Mesh{
public:
    Mesh(vector<int>Size);
    int GetDim(void);
    vector <int> GetSize(void);
    int GetNpoints(void);
private:
    int dim, Npnts;
    vector <int> size;
};

class Patch: public Mesh{
public:
    Patch(vector <int> Size, vector <double> bounds);
    vector<double> GetBound(void);
    vector<int> GetSteps(void);
    void ComputeCoords(const int i);
    vector<double> GetSpacing(void);
    vector < double> GetCoord(const int i);
private:
    vector <double> lengths;
    vector<vector <double>> coords;
    vector <double> spacing;
    vector <int> StencilSteps;
};

template <typename T> class DataMesh: public Mesh{
private:
    vector <T> field;
    double dt;
    int Nsteps;
public:
    DataMesh<T>(vector <int> Size, T f);
    DataMesh(vector <int> Size, const double Courant, const int Nstep);
    DataMesh<T>(vector <int> Size);
    ~DataMesh<T>(void);
    DataMesh<T> operator+(const DataMesh<T>& a);
    void operator += (const DataMesh<T>& b);
    void operator * (const T a);
    void SetElement (const int i, const T value);
    void EvolveField(void);
    void Print(void);
};

template <typename T> class ComputeRHS{
private:
    double dt;
    int Npnts;
    vector <int> StencilSteps;
    vector <int> Sizes;
    vector <T> derivative;
public:
    ComputeRHS<T>(vector <int> Size, const double Courant);
    ComputeRHS<T>(vector <int> Size);
    vector<T> CenteredDerivative1st (const vector <T>& field);
    vector<T> DownstreamDerivative1st (const vector<T>& field);
    vector<T> UpstreamDerivative1st (const vector<T>& field);
};



#endif //HW3_HW3_H
