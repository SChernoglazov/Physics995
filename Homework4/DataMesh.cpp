#include "hw.h"

template<typename T>
DataMesh<T>::DataMesh(){}

template<typename T>
DataMesh<T>::DataMesh(const vector <int> Size,const T f): Mesh(Size){
  int Npnt=GetNpoints();
  for (int i=0; i<Npnt; i++){
    field.push_back(f);
  }
  int step=1;
  for (int i=0; i<Size.size(); i++){
    StencilSteps.push_back(step);
    step*=Size[i];
  }
}

template<typename T>
DataMesh<T>::DataMesh(const vector <int> Size): Mesh(Size){
  int Npnt=GetNpoints();
  for (int i=0; i<Npnt; i++){
    field.push_back(0);
  }
  int step=1;
  for (int i=0; i<Size.size(); i++){
    StencilSteps.push_back(step);
    step*=(Size[i]);
  }
}

template<typename T>
DataMesh<T>::DataMesh(const int GhostZone, const vector <int> Size): Mesh(Size, GhostZone){
  int Npnt=GetNpoints();
  for (int i=0; i<Npnt; i++){
    field.push_back(0);
  }
  int step=1;
  for (int i=0; i<Size.size(); i++){
    StencilSteps.push_back(step);
    step*=(Size[i]+2*GhostZone);
  }
}

template <typename T>
void DataMesh<T>::SetGZMask(const int GhostZone, const vector <int> Size){ // the function is used to set GZ mask
  int Npnt=GetNpoints();
  int l, position;
  vector <int> reverse_coords;
  for (int i=0; i<Npnt; i++){
    l=i;
    for (int j=1; j<=Size.size();j++){
      position = l/StencilSteps[Size.size()-j];
      if (position < GhostZone || position >= (Size[Size.size()-j]+GhostZone))
        field[i]=1;
      l=l%StencilSteps[Size.size()-j];
    }
  }
}

template <typename T>
void DataMesh<T>::SetFaceGZ(const DataMesh<T>& GZ, const  int i){ //i is direction of the flux 0 -x, 1-y,2-z
  int Npnt=GetNpoints();
  bool flag;
  for (int j=0; j<Npnt; j++){
    flag=1;
    if (GZ.return_element(j)==0){
      flag=0;
    }
    if(j>=StencilSteps[i]){
      if(GZ.return_element(j-StencilSteps[i])==0){
	flag=0;
      }
    }
    field[j]=flag;
  }
}

template<typename T>
DataMesh<T>::~DataMesh(){}

template<typename T>
DataMesh<T> DataMesh<T>::operator+(const DataMesh<T>& a)
{
  vector<int> size=this->GetSize();
  if (a.field.size() != this->field.size()){
    cout << "a and b should have same size" << endl;
    exit(1);
  }

  DataMesh<T> c(size);
  if(typeid(T)== typeid(bool)) {
    for (int i = 0; i < a.field.size(); i++) {
      c.field[i] = a.field[i]*this->field[i];
    }
  }
  else{
    for (int i = 0; i < a.field.size(); i++) {
      c.field[i] = a.field[i] + this->field[i];
    }
  }
  return c;
}

template<typename T>
void DataMesh<T>::SetValue(const int i, const T& a){
  field[i]=a;
}
template <typename T>
void DataMesh<T>::operator += (const DataMesh<T>& b){
  if (b.field.size() != this->field.size()){
    cout << "a and b should have same size" << endl;
    exit(1);
  }
  if(typeid(T)== typeid(bool)){
    for (int i = 0; i < field.size(); i++) {
      field[i] = field[i] * b.field[i];
    }
  }
  else{
    for (int i=0; i<field.size(); i++){
      field[i]=field[i]+b.field[i];
    }
  }
}

template <typename T>
void DataMesh<T>::operator * (const T a){
  for (int i=0; i<field.size(); i++){
    field[i]=a*field[i];
  }
}

template <typename T>
T DataMesh<T>::return_element (const int i) const{
  return field[i];
}

template <typename T>
void DataMesh<T>::Print() const{
  for (int i=0; i<field.size(); i++){
    cout << "i=" << i << ", value=" << field[i] << endl;
  }
}
