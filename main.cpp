#include<iostream>
#include<cmath>
#include<vector>
#include "hw1.h"

using namespace std; 

double smooth(const double x) {
    double f;
    f = (x * x - 1) / (x + 2);
    return f;
}

double jump(const double x) {
    double f = 0;
    if (x > 0.15)
        f = 1;
    return f;
}

int main(){
    int n;
    int N;
    double dx, y, x, real_y, deviation;
    double *points;
    vector<double> resolution{0.8,0.4,0.2,0.15,0.1,0.075,0.05,0.04,0.025,0.0125,0.00625};
    cout << "if you want to see solution of problem 2 - press 2" << endl;
    cout << "problem 3 - press 3, and problem 4 - press 4" << endl;
    cin >> n;
    switch (n-1){
        case 1:
            N=5;
            dx=0.1;
            points = new double[N];
            for (int i = 0; i < N; i++) {
                points[i] = smooth(i * dx);
            }
            x = 2.5 * dx;
            y = interpolation(x, points, N, dx);
            cout << y << endl;
            break;
        case 2:
            N=5;
            points = new double[N];
            cout << "smooth function" << endl;
            for (int j = 0; j<resolution.size();j++) {
                dx = resolution[j];
                for (int i = 0; i < N; i++) {
                    points[i] = smooth(i * dx);
                }
                x = 2.5 * dx;
                y = interpolation(x, points, N, dx);
                real_y = smooth(x);
                deviation = fabs(y-real_y);
                cout << "resolution = " << dx << " error = " << deviation << endl;
            }
	    cout << "jump function" << endl;
	    for (int j = 0; j<resolution.size();j++) {
	      dx = resolution[j];
	      for (int i = 0; i < N; i++) {
       		points[i] = jump(i * dx);                                                                                                                                                                                                                               
	      }
	      x = 2.5 * dx;
	      y = interpolation(x, points, N, dx);
	      double real_y = jump(x);                                                                                                                                                                                                                                    
	      deviation = fabs(y-real_y);
	      cout << "resolution = " << dx << " error = " << deviation << endl;
            }
            break;
        case 3:
            N=5;
            dx=0.1;
            points = new double[N];
            for (int i = 0; i < N; i++) {
                points[i] = smooth(i * dx);
            }
            y = interpolation(points, N);
            cout << "Exact result for N=4 is -0.454651, and for N=5 is -0.41666" << endl;
            cout << "test for N=" << N << " result = " << y << endl;
	    N=4;
	    delete [] points;
            points = new double[N];
            for (int i = 0; i < N; i++) {
	      points[i] = smooth(i * dx);
            }
            y = interpolation(points, N);
            cout << "Exact result for N=4 is -0.454651, and for N=5 is -0.41666" << endl;
            cout << "test for N=" << N << " result = " << y << endl;
            break;
        default:
            break;
    }
    delete [] points;
    return 0;
}
