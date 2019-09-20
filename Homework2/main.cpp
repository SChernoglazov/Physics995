#include "hw2.h"
#include "DataMesh.cpp"

int main() {
    DataMesh<double> a({2,2,2},1), b({2,2,2},2);
    DataMesh<double> c=a+b;
/*    DataMesh<bool> a({2,2,2},0), b({2,2,2},0);
    DataMesh<bool> c=a+b;*/
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
    cout << "For the output I assume that the domain is three-dimensional" << endl;
    for (int i=0; i<d.GetNpoints(); i++)
        //cout << d.GetCoord(i) << endl;
        cout << "i=" << i << " (" << d.GetCoord(i)[0] << " " << d.GetCoord(i)[1] << " " << d.GetCoord(i)[2] << ")" << endl;
    return 0;
}