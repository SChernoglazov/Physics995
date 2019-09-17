
#ifndef HW2_HW2_H
#define HW2_HW2_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>

using namespace std;

class Mesh{
private:
    int dim, Npnts;
    vector <int> size;
public:
    Mesh(vector<int>Size);
    int GetDim(void);
    vector <int> GetSize(void);
    int GetNpoints(void);
};

template <typename T> class DataMesh: public Mesh{
private:
    vector <T> field;
public:
    DataMesh<T>(vector <int> Size, T f);
    DataMesh<T>(vector <int> Size);
    ~DataMesh<T>();
    DataMesh<T> operator+(const DataMesh<T>& a);
    void operator += (const DataMesh<T>& b);
    void operator * (const T a);
    void Print();
};

class Patch: public Mesh{
private:
    vector <double> lengths;
    vector<vector <double>> coords;
    vector <double> spacing;
    vector <int> StencilSteps;
public:
    Patch(vector <int> Size, vector <double> bounds);
    vector<double> GetBound(void);
    vector<int> GetSteps(void);
    void ComputeCoords(const int i);
    vector < double> GetCoord(const int i);
};

#endif //HW2_HW2_H
