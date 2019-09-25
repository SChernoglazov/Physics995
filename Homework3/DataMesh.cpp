#include "hw3.h"


template <typename T>
DataMesh<T>::DataMesh(vector <int> Size, T f): Mesh(Size){
    int Npnt=GetNpoints();
    for (int i=0; i<Npnt; i++){
        field.push_back(f);
    }
}

template <typename T>
DataMesh<T>::DataMesh(vector <int> Size): Mesh(Size){
    int Npnt=GetNpoints();
    for (int i=0; i<Npnt; i++){
        field.push_back(0);
    }
}

template <typename T>
DataMesh<T>::~DataMesh(){}

template <typename T>
DataMesh<T> DataMesh<T>::operator +(const DataMesh<T>& a)
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

template <typename T>
void DataMesh<T>::operator +=(const DataMesh<T>& b){
    if (b.field.size() != this->field.size()){
        cout << "a and b should have same size" << endl;
        exit(1);
    }
    if(typeid(T)== typeid(bool)){
        for (int i = 0; i < field.size(); i++) {
            this->field[i] = this->field[i] * b.field[i];
        }
    }
    else{
        for (int i=0; i<field.size(); i++){
            this->field[i]=this->field[i]+b.field[i];
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
void DataMesh<T>::SetElement(const int i, const T value){
    if (i>=field.size()){
        cout << "the index is more than the size of the field " << field.size() << endl;
        exit(1);
    } else {
        field[i]=value;
    }
}

template <typename T>
T DataMesh<T>::GetElement (const int i){
    return field[i];
}

template <typename T>
void DataMesh<T>::Print(){
    for (int i=0; i<field.size(); i++){
        cout << "i=" << i << ", value=" << field[i] << endl;
    }
}

