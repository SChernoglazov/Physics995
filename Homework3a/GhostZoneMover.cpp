#include "hw2.h"
#include "DataMesh.cpp"

GhostZoneMover::GhostZoneMover(vector<int> Size, int GhostZone){ 
  dim=Size.size();
  InternalSize=Size;
  Npnts=1;
  int steps=1;
  for (int i=0; i<dim;i++){
    Npnts*=(Size[i]+2*GhostZone);
    StencilSteps.push_back(steps);
    steps*=(Size[i]+2*GhostZone);
  }
}

void GhostZoneMover::PeriodicGZ(DataMesh<bool>& GZ, DataMesh<double>& field){
  double value;
  for (int i=0; i<Npnts; i++){
    if(GZ.return_element(i)==1){ //the point is inside the mask
      for (int j=0;j<dim; j++){
	if (i+InternalSize[j]*StencilSteps[j]<Npnts){
	  if (GZ.return_element(i+InternalSize[j]*StencilSteps[j])==0){//the point is outside of the GZ mask
	    value=field.return_element(i+InternalSize[j]*StencilSteps[j]);
	    field.SetValue(i, value);
	  }
	}
	if(i-InternalSize[j]*StencilSteps[j]>0){
	  if (GZ.return_element(i-InternalSize[j]*StencilSteps[j])==0){
	    value=field.return_element(i-InternalSize[j]*StencilSteps[j]);  
	    field.SetValue(i,value);
	  }
	}
      }
    }
  }
}

void GhostZoneMover::PeriodicGZ(DataMesh<bool>& GZ, vector<double>& field){
  double value;
  for (int i=0; i<Npnts; i++){
    if(GZ.return_element(i)==1){ //the point is inside the mask  
      for (int j=0;j<dim; j++){
        if (i+InternalSize[j]*StencilSteps[j]<Npnts){
          if (GZ.return_element(i+InternalSize[j]*StencilSteps[j])==0){//the point is outside of the GZ mask  
            value=field[i+InternalSize[j]*StencilSteps[j]];
            field[i]= value;
          }
        }
        if(i-InternalSize[j]*StencilSteps[j]>0){
          if (GZ.return_element(i-InternalSize[j]*StencilSteps[j])==0){
            value=field[i-InternalSize[j]*StencilSteps[j]];
            field[i]=value;
          }
        }
      }
    }
  }
}
