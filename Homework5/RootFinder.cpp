#include <iostream>
#include <cmath> 

using namespace std;

double f(double x){
  return cos(2/x)-2*sin(1/x)+1/x;
}

double Newtonian (double (*fun)(double x), double initial, double precision){
  double dx, df, d2f, error, x ,x1;
  dx = 0.01;
  x = initial;
  df = (fun(x+dx)-fun(x-dx))/dx;
  d2f = (fun(x-dx)+fun(x+dx) - 2*fun(x))/(2*dx);
  x1 = x - fun(x)/df;
  error = x1 - x;
  x=x1;
  while (abs(error) > precision){
    x1 = x - fun(x)/df;
    df = (fun(x+dx)-fun(x-dx))/dx;
    d2f = (fun(x-dx)+fun(x+dx) - 2*fun(x))/(2*dx);
    error = x1 - x;
    x=x1;
  }
  return x; 
}

double Bisection (double (*fun)(double x), double a, double b, double precision){
  if (fun(a)*fun(b)>0){
    cout << "the values of the function at two ends must be different" << endl;
    return 0;
  }
  double c, x, diff;
  x = (a+b)/2;
  diff = fabs(a-b);
  while (diff > precision/2){
    c = (a+b)/2;
    if (fun(a)*fun(c) <= 0){
      b = c;
    }else{
      a = c;
    }
    diff = fabs(a-b);
  }
  x = (a+b)/2;
  return x;
}

int main(){
  cout << "Newtonian " << Newtonian(f, 0.6, 0.000001) << endl;
  cout << "Bisectional " << Bisection(f, 0.1, 0.5, 0.0000001) << endl; 
  return 0;
}
