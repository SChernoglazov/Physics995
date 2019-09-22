#ifndef HW3_HW3_H
#define HW3_HW3_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <math.h>

#define PI 3.141592653589793
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
public:
    DataMesh<T>(vector <int> Size, T f);
    DataMesh<T>(vector <int> Size);
    ~DataMesh<T>(void);
    DataMesh<T> operator+(const DataMesh<T>& a);
    void operator += (const DataMesh<T>& b);
    void operator * (const T a);
    void SetElement (const int i, const T value);
    T GetElement (const int i);
    void Print(void);
};


template <typename T> class ComputeRHS{
private:
    int Npnts;
    vector <T> derivative;
    vector <int> StencilSteps;
public:
    ComputeRHS(vector<int> N);
    void UpstreamDerivative( DataMesh<T>& U, const double dt, DataMesh<T>& dtU);
    void DownstreamDerivative( DataMesh<T>& U, const double dt, DataMesh<T>& dtU);
    void CenteredDerivative( DataMesh<T>& U, const double dt, DataMesh<T>& dtU);
    void RungeKutta3 (DataMesh<T>& U, const double dt, DataMesh<T>& dtU);
    T ThirdDerivative(DataMesh<T>& U, const int i);
    T ThirdDerivative(vector<T>& U, const int i);
};



#endif //HW3_HW3_H
