#include "hw.h"

Mesh::Mesh(){}

Mesh::Mesh(const vector<int>Size){
    dim=Size.size();
    size=Size;
    Npnts=1;
    for (int i=0; i<Size.size();i++)
        Npnts *= Size[i];
}

Mesh::Mesh(const vector<int> Size, const int GhostZone){
  dim=Size.size();
  Npnts=1;
  for (int i=0; i<dim; i++){
    Npnts *=(Size[i]+2*GhostZone);
    size.push_back(Size[i]+2*GhostZone);
  }
}

int Mesh::GetDim() const{
    return dim;
}
vector <int> Mesh::GetSize() const{
    return size;
}
int Mesh::GetNpoints() const{
    return Npnts;
}
