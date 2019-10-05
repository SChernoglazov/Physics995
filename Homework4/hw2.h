
#ifndef HW2_HW2_H
#define HW2_HW2_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <cmath>
#define PI 3.141592653589793

using namespace std;

class Mesh{
private:
    int dim, Npnts;
    vector <int> size;
public:
    Mesh();
    Mesh(const vector<int>Size);
    Mesh(const vector<int>Size,const int GhostZone);
    int GetDim(void)const;
    vector <int> GetSize(void)const;
    int GetNpoints(void)const;
};

template <typename T> class DataMesh: public Mesh{
private:
    vector <T> field;
    vector <int> StencilSteps;
public:
    DataMesh<T>();
    DataMesh<T>(const vector <int> Size, const T f);
    DataMesh<T>(const vector <int> Size);
    DataMesh<T>(const int GhostZone, const vector <int> Size);
    void SetGZMask(const int GhostZone, const vector <int> Size);
    void SetFaceGZ(const DataMesh<T>& GZ, const int i);
    ~DataMesh<T>();
    DataMesh<T> operator+(const DataMesh<T>& a);
    void operator += (const DataMesh<T>& b);
    void SetValue(const int i, const T& a);
    void operator * (const T a);
    T return_element(const int i)const;
    void Print()const;
};

class Patch: public Mesh{
private:
    vector <double> lengths;
    vector <DataMesh<double>> coords;
    vector <double> spacing;
    vector <int> StencilSteps;
public:
    Patch(const vector <int> Size, const  vector <double> bounds);
    Patch(const vector <int> Size, const  vector <double> bounds, const int GhostZone);
    vector<double> GetBound(void) const;
    vector<int> GetSteps(void) const;
    void ComputeCoords(const int i);
    void ComputeCoords(const int i, const int GhostZone);
    vector<double> GetCoord(const int i);
};

class GhostZoneMover{
private:
  int dim, Npnts, GhostZones;
  vector<int> InternalSize;
  vector<int> StencilSteps;
  vector<vector<int>> map;
public:
  GhostZoneMover(const vector<int> Size, const int GhostZone);
  void GeneratePeriodicGZ(const DataMesh<bool>& GZ,const DataMesh<double>& field);
  void ApplyBCs(DataMesh<double>& field);
  void ApplyBCs(vector<double>& field);
};

template <typename T> class ComputeRHS{
 private:
  int Npnts, GhostZones;
  vector <int> StencilSteps;
  vector<int> Sizes;
 public:
  ComputeRHS();
  void FillComputeRHS(const vector<int> N, const int GhostZone);
  T ThirdDerivative(const DataMesh<T>& U, const int i, const DataMesh<bool>& GZ);
  T ThirdDerivative(const vector<T>& U, const int i, const DataMesh<bool>& GZ);
  T BurgersFlux(const T u) const;
  T AdvectionFlux(const T u) const;
};

template <typename T> class  dtU{
 private:
  int Npnts, GhostZones, dim;
  vector <int> StencilSteps;
  vector<int> Sizes;
  vector<vector<T>> FluxLR;
  vector<vector<T>> Flux;
  vector<vector<T>> ULR;
  ComputeRHS<T> RHS;
 public:
  dtU(const vector <int> Size, const int GhostZone);
  void RungeKutta3(const DataMesh<T>& U, const double dt, DataMesh<T>& dUt, const DataMesh<bool>& GZ, GhostZoneMover& GZM);
  void RungeKutta3shock(const DataMesh<T>& U, const double dt, DataMesh<T>& dUt, const DataMesh<bool>& GZ, GhostZoneMover& GZM, const vector<DataMesh<bool>>& FaceGZ);
  void FluxLRReconstructor(const vector<DataMesh<bool>>& FaceGZ, const DataMesh<T>& U);
  void FluxLRReconstructor(const vector<DataMesh<bool>>& FaceGZ, const vector<T>& U);
  void FluxReconstructor(const vector<DataMesh<bool>>& FaceGZ);
  void divFlux(const DataMesh<bool>& GZ,vector<T>& helper);
  int sign(const T a) const;
  T MinMod(const T a, const T b, const T c) const;
};

#endif //HW2_HW2_H
