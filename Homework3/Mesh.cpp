#include "hw3.h"

Mesh::Mesh(vector<int>Size){
    dim=Size.size();
    size=Size;
    Npnts=1;
    for (int i=0; i<Size.size();i++)
        Npnts *= Size[i];
}

int Mesh::GetDim(){
    return dim;
}

vector <int> Mesh::GetSize(){
    return size;
}
int Mesh::GetNpoints(){
    return Npnts;
}
