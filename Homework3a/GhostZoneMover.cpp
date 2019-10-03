#include "hw2.h"
#include "DataMesh.cpp"

GhostZoneMover::GhostZoneMover(vector<int> Size, int GhostZone){ 
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

void GhostZoneMover::PeriodicGZ(DataMesh<bool>& GZ, DataMesh<double>& field){ 
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
    	  value=field.return_element(i+InternalSize[j]*StencilSteps[j]);
	  field.SetValue(i, value);
	}
	if(coords[j]>=InternalSize[j]+GhostZones){
	  value=field.return_element(i-InternalSize[j]*StencilSteps[j]);
          field.SetValue(i, value);
	}
      }
    }
  }
}

void GhostZoneMover::PeriodicGZ(DataMesh<bool>& GZ, vector<double>& field){
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
          value=field[i+InternalSize[j]*StencilSteps[j]];
          field[i] = value;
        }
        if(coords[j]>=InternalSize[j]+GhostZones){
          value=field[i-InternalSize[j]*StencilSteps[j]];
          field[i] = value;
        }
      }
    }
  }
}

