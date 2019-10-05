#include "hw2.h"
#include "DataMesh.cpp"

GhostZoneMover::GhostZoneMover(const vector<int> Size, const int GhostZone){ 
  dim=Size.size();
  InternalSize=Size;
  Npnts=1;
  GhostZones=GhostZone;
  int steps=1;
  for (int i=0; i<dim;i++){
    Npnts*=(Size[i]+2*GhostZone);
    StencilSteps.push_back(steps);
    steps*=(Size[i]+2*GhostZone);
  }
}

void GhostZoneMover::GeneratePeriodicGZ(const DataMesh<bool>& GZ, const DataMesh<double>& field){ 
  vector <int> coords(dim);
  int l, step;
  double value=0;
  for (int i=0; i<Npnts; i++){
    if (GZ.return_element(i)==1){ // if the points is inside the GZ
      l=i;
      for (int j=1; j<=dim; j++){ //recover the coordinates
	step = l/StencilSteps[dim-j];
        coords[j-1]=step;
        l=l%StencilSteps[dim-j];
      }
      reverse(coords.begin(),coords.end());
      for (int j=0; j<dim; j++){
	if(coords[j]<GhostZones){
	  map.push_back({i,i+InternalSize[j]*StencilSteps[j]});
	  continue;
	}
	if(coords[j]>=InternalSize[j]+GhostZones){
	  map.push_back({i,i-InternalSize[j]*StencilSteps[j]});
	}
      }
    }
  }
}

void GhostZoneMover::ApplyBCs(DataMesh<double>& field){
  double value;
  for (int i=0; i<map.size(); i++){
    value=field.return_element(map[i][1]);
    field.SetValue(map[i][0],value);
  }
}

void GhostZoneMover::ApplyBCs(vector<double>& field){
  double value;
  for (int i=0; i<map.size(); i++){
    value=field[map[i][1]];
    field[map[i][0]]=value;
  }
}
