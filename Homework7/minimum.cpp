#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

vector<double> p;
vector<double> n;

double f(double x)
{
  double alpha=-0.143*sin(1.75*(x+1.73));
  double beta=-0.18*sin(2.96*(x+4.98));
  double delta=0.012*sin(6.23*(x+3.17));
  double gamma=0.088*sin(8.07*(x+4.63));
  double val=alpha+beta+delta+gamma;
  return val;
}

int sign(double x){
  return (x>=0 ? 1 : -1);
}

double max (double a, double b){
  return (a>b ? a : b);
}

vector<double> bracketing_minimum(double (*fun)(double x), double initial){
  double a, b, c, dx=0.1, aux;
  double step=1.618; //golden section
  int count = 0;
  a=initial;
  b=initial+dx;
  if(fun(b)>fun(initial)){
    b=initial;
    a=initial+dx;
    dx=-dx;
  }
  dx = dx*step;
  c=b+dx;
  //cout << " a =" << a <<" b=" << b << " c=" << c <<  endl;
  while((fun(c)-fun(b))*(fun(b)-fun(a))>=0){
    dx = dx*step;
    aux=c+dx;
    a=b;
    b=c;
    c=aux;
    //    cout << " a =" << a <<" b=" << b << " c=" << c <<  endl;
    count++;
    if(count>=100){
      cout << "too many itterations for bracketing" << endl;
      break;
      }
  }
  if(count>=100)
    return {0,0};
  else
    return {min(a,c),b,max(a,c)};
}


double golden_section (double (*fun)(double x),double a, double b){
  double tol = 3.e-8;; //maximum tolerance level for golden section method
  double R=0.618003399, w=1-R;
  double xleft, xright, xmid, xnew;
  int count=0;
  if (a>b){//to make sure c>a
    xright=a;
    xleft=b;
  }else{
    xright=b;
    xleft=a;
  }
  xmid=xleft+(xright-xleft)*w;
  double err=xright-xleft;
  while(err>tol){
    if((xmid-xleft) > (xright-xmid)){ // if left subinterval is larger
      xnew=xleft+(xmid-xleft)*w;
      if (fun(xnew)>fun(xmid)){
	xleft=xnew;
      }else{
	xright=xmid;
	xmid=xnew;
      }
    }else{ //if right one is larger
      xnew=xmid+(xright-xmid)*w;
      if (fun(xnew)>fun(xmid)){
	xright=xnew;
      }else{
	xleft=xmid;
	xmid=xnew;
	}
    }
    err=xright-xleft;
    count++;
    if(count>=100){//maximum number of iterations 
      return a;
    }
  }
  return (xleft+xright)/2;
}
  
double Brendt(double (*fun)(double x),double ax, double bx, double cx){
  double xleft, xright, x, u, w, v;
  double r,q,p, fu, fx, fv, fw;
  double e=0, etemp, tol=3e-8, tol2, d;
  double R=0.618003399, G=1-R;
  int count=0;
  double xm; //midpoint
  xleft=(ax<cx ? ax : cx);
  xright=(ax>cx? ax : cx);
  // the criterium of convergence is all new points are in between a and b
  x=w=v=u=bx;
  fw=fv=fx=fun(x);
  xm=(xleft+xright)/2;
  while (abs(x-xm)>(tol2-0.5*(xright-xleft))){
    xm=(xleft+xright)/2;
    tol2=2*abs(x)*tol;
    if (count==100){
      cout << " too many itterations" << endl;
      return xleft;
    }
    count++;
    if(abs(e)>tol*abs(x)){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;//nominator
      q=2*(q-r);//denominator
      if (q>0)
	p=-p;
	q=abs(q);
	etemp=e;
	e=d;
	if(abs(p) >= abs(0.5*q*etemp) || p <= q*(xleft-x) || p >=q*(xright-x)){
	  // first condition - step is larger than one before last
	  // 2nd & 3rd condition - parabolic fit puts the new point outside of the bracket
	  // apply golden section to the larger part
	  d = G*(x>xm ? xleft-x : xright-x);
	  e = (x>xm ? xleft-x : xright-x);
	}else{
	  d = p/q;
	  u=x+d;//parabolic fit
	  if ((u-xleft) < tol2 || (xright-u) <tol2){
	    //new point is very close to the boundary of the domain
	    // and it is inside the brackets
	    //useless to move less than tolerance of the method 
	    d= abs(x)*tol;
	  }
	}
    }else{
      d = G*(x>xm ? xleft-x : xright-x);
      e = (x>xm ? xleft-x : xright-x);
    }
    //useless to move less than tolerance of the method
    u=(abs(d)>=abs(x)*tol ? x+d : x+sign(d)*abs(x)*tol);
    fu=fun(u);
    if (fu<fx){ //are we going down towards minimum?
      // what interval we get after the itteration
      if (u >= x)
	xleft=x;
      else 
	xright=x;
      v=w;
      w=x;
      x=u;
      fv=fw;
      fw=fx;
      fx=fu;
    }else{
      if (u<x)
	xleft=u;
      else 
	xright=u;
      if(fu<=fw || w == x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }else if(fu <=fv || v==x || v==w){
	v=u;
	fv=fu;
      }
    }
  }
  return x;
}

double f2(vector<double>& arg){
  double res=1;
  for (int i=0; i<arg.size(); i++)
    res*=pow(arg[i],4);
  return res;
}

double falong (double lambda){
  vector<double> x;
  for (int i=0; i<p.size(); i++)
    x.push_back(p[i]+lambda*n[i]);
  return f2(x);
}

void find_min_along(double initial){
  vector<double> bracket = bracketing_minimum(falong, initial);
  cout << "2d bracketing for parametrization parameter is (" 
       << bracket[0] <<  ", " << bracket[1] << ", " << bracket[2] <<")"<< endl;
  cout << "2d golden section minimum is for parametrization along a line " << golden_section(falong,bracket[0],bracket[2]) << endl;
  cout << "2d Brendt's method minimum is for parametrization along a line " << Brendt(falong,bracket[0],bracket[1],bracket[2]) << endl;
}


vector<double> steepest_descent (){
  vector <double> oldp=p;
  vector <double> oldn=n;
  vector <double> gradient (p.size());
  vector <double> new_point, bracket;
  vector <double> positive_shift, negative_shift;
  double dx=0.001, tol=3e-8, diff=2*tol, lambda;
  int count=0;
  bool flag;
  new_point=p;
  while (diff>tol){
    diff=0;
    //computation of the gradient
    for (int i=0; i<p.size(); i++){
      positive_shift=negative_shift=new_point;
      positive_shift[i] += dx;
      negative_shift[i] -= dx;
      //minus gradient
      gradient[i]= (f2(negative_shift)-f2(positive_shift))/(2*dx);
    }
    flag=0;
    // if the gradient is effectively 0, we reached some minimum... 
    for (int i=0; i<p.size(); i++){
      if(fabs(gradient[i])>1e-16){
	flag=1;
      }
    }
    if (!flag){
      break;
    }
    n=gradient;
    p=new_point;
    // finding the mimumum along the line
    bracket = bracketing_minimum(falong, 0);
    lambda = Brendt(falong,bracket[0],bracket[1],bracket[2]);
    for (int i=0; i<p.size(); i++){
      new_point[i]+= lambda*n[i]; 
      // how large is the step? If it small, we've reached minimum
      diff += pow(lambda*n[i],2);
    }
    diff=sqrt(diff)/p.size(); 
    count++;
    if (count == 100){
      cout << "too many itterations of the steepest descent method " << endl;
      break;
    }
  }
  //recovery old p and n
  p=oldp;
  n=oldn;
  if (count < 100)
    return new_point;
  else
    return oldp;
}

vector<double> grad(vector<double>& x){
  vector <double> positive_shift, negative_shift;
  vector <double> gradient(x.size());
  double dx=0.001;
  for (int i=0; i<x.size(); i++){
    positive_shift=negative_shift=x;
    positive_shift[i] += dx;
    negative_shift[i] -= dx;
    gradient[i]= (f2(positive_shift)-f2(negative_shift))/(2*dx);
  }
  return gradient;
}

vector<double> conjugate_gradient (){
  double eps=1e-16, tol=1e-8, ftol=3e-8;
  double gg, dgg, test, den, temp, gam, fret, lambda;
  int N=p.size(), i;
  vector<double> pp=p;
  vector<double> nn=n;
  vector<double> new_point;
  vector<double> bracket;
  vector<double> g(N), h(N);
  double fp=f2(p);
  n=grad(p);
  for(int j=0; j<N; j++){
    g[j]=-n[j];
    n[j]=h[j]=g[j];
  }
  for (i=0; i<200; i++){
    bracket = bracketing_minimum(falong, 0);
    lambda = Brendt(falong,bracket[0],bracket[1],bracket[2]);
    for (int j=0; j<N; j++){
      p[j] +=lambda*n[j];
    }
    fret=f2(p);
    if(2*fabs(fret-fp) <=ftol*(fabs(fret)+fabs(fp)+eps)){
      new_point=p;
      break;
    }
    fp=fret;
    n=grad(p);
    test=0;
    den=max(fp, 1.0);
    for (int j=0; j<N;j++){
      temp=fabs(n[j])*max(fabs(p[j]),1.0)/den;
      if (temp>test)
	test=temp;
    }
    if (test<tol){
      new_point=p;
      break;
    }
    dgg=gg=0;
    for (int j=0; j<N; j++){
      gg += g[j]*g[j];
      dgg +=(n[j]+g[j])*n[j];
    }
    if (gg<1e-31){ // gradient is zero
      new_point=p;
      break;
    }
    gam=dgg/gg;
    for (int j=0; j<N; j++){
      g[j] = -n[j];
      n[j]=h[j]=g[j]+gam*h[j];
    }
  }
  if(i>=199){
    cout << "too many itterations for conjugate gradient"  << endl;
    return pp;
  }else{
    //recovery old values
    p=pp;
    n=nn;
    return new_point; 
  }
}

int main(){
  vector<double> bracket;
  p = {1,-1};
  n = {0,1};
  bracket=bracketing_minimum(f, 0.3);
  if (bracket[0]==bracket[1])
    cout<<"bracketing failed" << endl;
  else
    cout << "bracketing is (" << bracket[0] <<  ", " << bracket[1] << ", " << bracket[2] <<")"<< endl;
  cout << "golden section minimum is " << golden_section(f,bracket[0],bracket[2]) << endl;
  cout << "Brendt's method" << Brendt(f,bracket[0],bracket[1],bracket[2]) << endl;
  find_min_along(-3);
  cout << "steepest descent " << steepest_descent()[0] <<" " <<  steepest_descent()[1] << endl;
  cout << "conjugate gradient method point " << conjugate_gradient()[0] << conjugate_gradient()[1] << endl;
  return 0;
}
