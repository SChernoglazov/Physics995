
#ifndef HW2_HW2_H
#define HW2_HW2_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#define PI 3.141592653589793

using namespace std;

class Mesh{
private:
    int dim, Npnts;
    vector <int> size;
public:
    Mesh(vector<int>Size);
    Mesh(vector<int>Size, int GhostZone);
    int GetDim(void);
    vector <int> GetSize(void);
    int GetNpoints(void);
};

template <typename T> class DataMesh: public Mesh{
private:
    vector <T> field;
    vector <int> StencilSteps;
public:
    DataMesh<T>(vector <int> Size, T f);
    DataMesh<T>(vector <int> Size);
    DataMesh<T>(int GhostZone, vector <int> Size);
    void SetGZMask(int GhostZone, vector <int> Size);
    ~DataMesh<T>();
    DataMesh<T> operator+(const DataMesh<T>& a);
    void operator += (const DataMesh<T>& b);
    void SetValue(const int i, const T& a);
    void operator * (const T a);
    T return_element(const int i);
    void Print();
};

class Patch: public Mesh{
private:
    vector <double> lengths;
    vector <DataMesh<double>> coords;
    vector <double> spacing;
    vector <int> StencilSteps;
public:
    Patch(vector <int> Size, vector <double> bounds);
    Patch(vector <int> Size, vector <double> bounds,int GhostZone);
    vector<double> GetBound(void);
    vector<int> GetSteps(void);
    void ComputeCoords(const int i);
    void ComputeCoords(const int i, const int GhostZone);
    vector<double> GetCoord(const int i);
};

class GhostZoneMover{
private:
  int dim, Npnts;
  vector<int> InternalSize;
  vector<int> StencilSteps;
public:
  GhostZoneMover(vector<int> Size, int GhostZone);
  void PeriodicGZ(DataMesh<bool>& GZ, DataMesh<double>& field);
  void PeriodicGZ(DataMesh<bool>& GZ, vector<double>& field);
};

template <typename T> class ComputeRHS{
 private:
  int Npnts, GhostZones;
  vector <int> StencilSteps;
  vector<int> Sizes;
 public:
  ComputeRHS();
  void FillComputeRHS(vector<int> N, int GhostZone);
  T ThirdDerivative(DataMesh<T>& U, const int i, DataMesh<bool>& GZ);
  T ThirdDerivative(vector<T>& U, const int i, DataMesh<bool>& GZ);
};

template <typename T> class  dtU{
private:
  int Npnts, GhostZones;
  vector <int> StencilSteps;
  vector<int> Sizes;
  ComputeRHS<T> RHS;
public:
  dtU(vector <int> Size, int GhostZone);
  void RungeKutta3(DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM);
  void UpstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM);
  void DownstreamDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM);
  void CenteredDerivative(DataMesh<T>& U, const double dt, DataMesh<T>& dUt, DataMesh<bool>& GZ, GhostZoneMover& GZM);
};

#endif //HW2_HW2_H
