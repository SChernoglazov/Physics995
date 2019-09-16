#include "hw3.h"

Patch::Patch(vector <int> Size, vector <double> bounds): Mesh(Size){
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
vector<double> Patch::GetBound(){
    return lengths;
}
vector<int> Patch::GetSteps(){
    return StencilSteps;
}

vector<double> Patch::GetSpacing(){
    return spacing;
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
    coords.push_back(reverse_coords);
}
vector < double> Patch::GetCoord(const int i){
    return coords[i];
}

