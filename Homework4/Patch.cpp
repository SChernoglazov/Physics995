#include "hw2.h"
#include "DataMesh.cpp"

Patch::Patch(const vector <int> Size, const vector <double> bounds): Mesh(Size){
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

Patch::Patch(const vector <int> Size, const vector <double> bounds, const int GhostZone): Mesh(Size, GhostZone){
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
      step*=(Size[k]+2*GhostZone);
    StencilSteps.push_back(step);
  }
  cout << StencilSteps[0] << " " << StencilSteps[1] << " " << StencilSteps[2] << endl;
  int Npnts=GetNpoints();
  for (int i=0; i<Npnts; i++){
    ComputeCoords(i, GhostZone);
  }
}

vector<double> Patch::GetBound() const{
    return lengths;
}
vector<int> Patch::GetSteps() const{
    return StencilSteps;
}

void Patch::ComputeCoords(const int i){ // number of the point
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
    int dim=reverse_coords.size();
    DataMesh<double> coordinate({dim});
    for (int i=0; i<dim; i++){
        coordinate.SetValue(i,reverse_coords[i]);
    }
    coords.push_back(coordinate);
}

void Patch::ComputeCoords(const int i, const int GhostZone){  
  int l=i, step,length;
  double local;
  vector <double> reverse_coords;
  vector<int> sizes = GetSize();
  length=sizes.size();
  for (int j=1; j<=length;j++){
    step = l/StencilSteps[sizes.size()-j];
    local=(step-GhostZone)*spacing[sizes.size()-j];
    reverse_coords.push_back(local);
    l=l%StencilSteps[sizes.size()-j];
  }
  reverse(reverse_coords.begin(),reverse_coords.end());
  int dim=reverse_coords.size();
  DataMesh<double> coordinate({dim});
  for (int i=0; i<dim; i++){
    coordinate.SetValue(i,reverse_coords[i]);
  }
  coords.push_back(coordinate);
}

vector < double> Patch::GetCoord(const int i){
    int dim = StencilSteps.size();
    vector<double> result;
    for (int j=0; j<dim; j++){
        result.push_back(coords[i].return_element(j));
    }
    return result;
}
