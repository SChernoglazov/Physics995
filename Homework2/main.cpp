#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Mesh{
private:
    int dim, Npnts;
    vector <int> size;
public:
    Mesh(vector<int>Size){
        dim=Size.size();
        size=Size;
        Npnts=1;
        for (int i=0; i<Size.size();i++)
            Npnts *= Size[i];
    }
    int GetDim(){
        return dim;
    }
    vector <int> GetSize(){
        return size;
    }
    int GetNpoints(){
        return Npnts;
    }
};

template <typename T> class DataMesh: public Mesh{
private:
    vector <T> field;
public:
    DataMesh<T>(vector <int> Size, T f): Mesh(Size){
        int Npnt=GetNpoints();
        for (int i=0; i<Npnt; i++){
            field.push_back(f);
        }
    }
    DataMesh<T>(vector <int> Size): Mesh(Size){
        int Npnt=GetNpoints();
        for (int i=0; i<Npnt; i++){
            field.push_back(0);
        }
    }
    ~DataMesh<T>(){}

    DataMesh<T> operator+(const DataMesh<int>& a)
    {
        vector<int> size=this->GetSize();
        if (a.field.size() != this->field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }

        DataMesh<T> c(size);
        for (int i=0; i<a.field.size(); i++){
            c.field[i]=a.field[i]+this->field[i];
        }
        return c;
    }

    DataMesh<T> operator+(const DataMesh<double>& a)
    {
        vector<int> size=this->GetSize();
        if (a.field.size() != this->field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }

        DataMesh<T> c(size);
        for (int i=0; i<a.field.size(); i++){
            c.field[i]=a.field[i]+this->field[i];
        }
        return c;
    }

    DataMesh<T> operator+(const DataMesh<bool>& a)
    {
        vector<int> size=this->GetSize();
        if (a.field.size() != this->field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }

        DataMesh<T> c(size);
        for (int i=0; i<a.field.size(); i++){
            c.field[i]=a.field[i]+this->field[i];
        }
        return c;
    }

    void operator += (const DataMesh<int>& b){
        if (field.size() != b.field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }
        for (int i=0; i<field.size(); i++){
            field[i]=field[i] +b.field[i];
        }
    }

    void operator += (const DataMesh<double>& b){
        if (field.size() != b.field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }
        for (int i=0; i<field.size(); i++){
            field[i]=field[i] +b.field[i];
        }
    }

    void operator += (const DataMesh<bool>& b){
        if (field.size() != b.field.size()){
            cout << "a and b should have same size" << endl;
            exit(1);
        }
        for (int i=0; i<field.size(); i++){
            field[i]=field[i] +b.field[i];
        }
    }

    void operator * (const double a){
        for (int i=0; i<field.size(); i++){
            field[i]=a*field[i];
        }
    }
    void Print(){
        for (int i=0; i<field.size(); i++){
            cout << "i=" << i << ", value=" << field[i] << endl;
        }
    }
};

class Patch: public Mesh{
private:
    vector <double> lengths;
    vector<vector <double>> coords;
    vector <double> spacing;
    vector <int> StencilSteps;
public:
    Patch(vector <int> Size, vector <double> bounds): Mesh(Size){
        lengths = bounds;
        int step;
        if ((bounds.size()-2*Size.size())!=0){
            cout << "size of vector of boundaries should be" << endl;
            cout << "twice larger than the dimension of the grid." << endl;
            exit(1);
        }
        StencilSteps.push_back(1);
        for (int i =0; i<Size.size(); i++){
            spacing.push_back((lengths[i*2+1]-lengths[i*2])/(Size[i]-1));
            step = 1;
            for (int k=0; k<=i; k++)
                step*=Size[k];
            StencilSteps.push_back(step);
        }
        int Npnts=GetNpoints();
        for (int i=0; i<Npnts; i++){
            ComputeCoords(i);
        }
    }
    vector<double> GetBound(){
        return lengths;
    }
    vector<int> GetSteps(){
        return StencilSteps;
    }
    void ComputeCoords(const int i){ // number of the point
        int l=i, step,length;
        double local;
        vector <double> reverse_coords;
        vector<int> sizes = GetSize();
        length=sizes.size();
        for (int j=1; j<=length;j++){
            step = l/StencilSteps[sizes.size()-j];
            local=step*spacing[sizes.size()-j];
            reverse_coords.push_back(local);
            l=l%StencilSteps[sizes.size()-j];
        }
        reverse(reverse_coords.begin(),reverse_coords.end());
        coords.push_back(reverse_coords);
    }
    vector < double> GetCoord(const int i){
        return coords[i];
    }
};

int main() {
    DataMesh<double> a({2,2,2},1), b({2,2,2},2);
    DataMesh<double> c=a+b;
    cout << "test of +: should be 3 everywhere" << endl;
    c.Print();
    cout << "test of +=: should be 3 everywhere" << endl;
    a+=b;
    a.Print();
    cout << "multiplication by number 6.7. As a=3" << endl;
    cout << "the result should be either 20 for int or 20.1 for double" << endl;
    a*6.7;
    a.Print();
    Patch d({4,6,2},{0,3,0,5,0,2});
    cout << "list of coordinates" << endl;
    for (int i=0; i<d.GetNpoints(); i++)
        cout << "i=" << i << " (" << d.GetCoord(i)[0] << " " << d.GetCoord(i)[1] << " " << d.GetCoord(i)[2] << ")" << endl;
    return 0;
}