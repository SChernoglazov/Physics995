#include "hw3.h"


template <typename T>
DataMesh<T>::DataMesh(vector <int> Size, T f): Mesh(Size){
    int Npnt=GetNpoints();
    for (int i=0; i<Npnt; i++){
        field.push_back(f);
    }
    dt=0.0025;
    Nsteps=200;
}

template <typename T>
DataMesh<T>::DataMesh(vector <int> Size): Mesh(Size){
    int Npnt=GetNpoints();
    for (int i=0; i<Npnt; i++){
        field.push_back(0);
    }
    dt=0.0025;
    Nsteps=200;
}

template <typename T>
DataMesh<T>::DataMesh(vector <int> Size, const double Courant, const int Nstep): Mesh(Size){
    int Npnt=GetNpoints();
    for (int i=0; i<Npnt; i++){
        field.push_back(0);
    }
    dt=Courant;
    Nsteps=Nstep;
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
void DataMesh<T>::SetElement(const int i, const T value){
    if (i>=field.size()){
        cout << "the index is more than the size of the field " << field.size() << endl;
        exit(1);
    } else {
        field[i]=value;
    }
}

template <typename T>
void DataMesh<T>::Print(){
    for (int i=0; i<field.size(); i++){
        cout << "i=" << i << ", value=" << field[i] << endl;
    }
}

template <typename T>
void DataMesh<T>::EvolveField(){
    int Npnt=GetNpoints();
    vector<int> size=this->GetSize();
    ComputeRHS<T> der(size, dt);
    for (int j=0; j<=Nsteps; j++) {
        if (j%40000 ==0) {
            cout << "N = "<< j << endl;
            cout << "{";
            for (int i = 0; i < Npnt-1; i++)
                cout << " "<<field[i] <<",";
            cout << " "<<field[Npnt-1];
            cout << "}" << endl;
        };
        vector<T> derivative = der.UpstreamDerivative1st(field);
        for (int i = 0; i < Npnt; i++) {
            field[i] += derivative[i];
        }
    }
}
