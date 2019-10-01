#include "hw2.h"

Mesh::Mesh(vector<int>Size){
    dim=Size.size();
    size=Size;
    Npnts=1;
    for (int i=0; i<Size.size();i++)
        Npnts *= Size[i];
}

Mesh::Mesh(vector<int> Size, int GhostZone){
  dim=Size.size();
  Npnts=1;
  for (int i=0; i<dim; i++){
    Npnts *=(Size[i]+2*GhostZone);
    size.push_back(Size[i]+2*GhostZone);
  }
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
