#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;

double Pres(double T, double Rho0){
  return 100*Rho0*Rho0+Rho0*T;
}

double Enth(double T, double Rho0){
  return 1+2*Pres(T, Rho0)/Rho0;
}

double RhoStar(double Detg, double W, double Rho0){
  return Detg*W*Rho0;
}

double Tau (double Detg, double W, double Rho0, double T){
  return RhoStar(Detg,W,Rho0)*(Enth(T,Rho0)*W-1)-Detg*Pres(T,Rho0);
}

vector<double> Sk (double Rho0, double T, double W, double Detg,vector<double>u_k){
  vector<double> result(u_k.size());
  for (int i=0; i<u_k.size(); i++)
    result[i]=RhoStar(Detg, W, Rho0)*Enth(T,Rho0)*u_k[i];
  return result;
}

class Inversion{
private:
  double S, RhoStar, Tau, Detg=1;
  vector <double> S_k;
  double precision=3*1e-16;
  double Wfound, Tfound;
  bool polynomial;
public:
  Inversion(double RhoSt, double tau, vector<double> Sk, bool poly){
    RhoStar=RhoSt;
    Tau=tau;
    S=0;
    for (int i=0; i<Sk.size(); i++){
      S_k.push_back(Sk[i]);
      S+=S_k[i]*S_k[i];
    }
    S = sqrt(S);
    polynomial=poly;
  }
  
  double Pres(double T, double Rho0){
    return 100*Rho0*Rho0+Rho0*T;
  }
  
  double Enth(double T, double Rho0){
    return 1+2*Pres(T, Rho0)/Rho0;
  }

  double Equation1 (double W, double T){
    double Rho0=RhoStar/(W*Detg);
    return W*W-1-S*S/pow(RhoStar*Enth(T,Rho0),2.0);
  }

  double Equation2 (double W, double T){
    double Rho0=RhoStar/(Detg*W);
    return Enth(T,Rho0)*W-1-Detg*Pres(T,Rho0)/RhoStar-Tau/RhoStar;
  }
  
  void findWT(double Wg, double Tg){ //Wg and Tg are guess for W and T
    double dW=1e-6, dT=1e-6;
    double f1W, f1T, f2W, f2T, f1, f2;// f1W = d f1(W,T)/d W
    double errorW, errorT, W, T;
    int count=0;
    f1W = (Equation1(Wg+dW, Tg)-Equation1(Wg-dW, Tg))/(2*dW);
    f2W = (Equation2(Wg+dW, Tg)-Equation2(Wg-dW, Tg))/(2*dW);
    f1T = (Equation1(Wg, Tg+dT)-Equation1(Wg, Tg-dT))/(2*dT);
    f2T = (Equation1(Wg, Tg+dT)-Equation1(Wg, Tg-dT))/(2*dT);
    f1 = Equation1(Wg, Tg);
    f2 = Equation2(Wg, Tg);
    errorW = (f2*f1T-f1*f2T)/(f1W*f2T-f2W*f1T);
    errorT = (f1*f2W-f2*f1W)/(f1W*f2T-f2W*f1T);
    Wg = Wg + errorW;
    Tg = Tg + errorT;
    while ((fabs(errorW)>precision)&&(fabs(errorT)>precision)){
      f1W = (Equation1(Wg+dW, Tg)-Equation1(Wg-dW, Tg))/(2*dW);
      f2W = (Equation2(Wg+dW, Tg)-Equation2(Wg-dW, Tg))/(2*dW);
      f1T = (Equation1(Wg, Tg+dT)-Equation1(Wg, Tg-dT))/(2*dT);
      f2T = (Equation1(Wg, Tg+dT)-Equation1(Wg, Tg-dT))/(2*dT);
      f1 = Equation1(Wg, Tg);
      f2 = Equation2(Wg, Tg);
      errorW = (f2*f1T-f1*f2T)/(f1W*f2T-f2W*f1T);
      errorT = (f1*f2W-f2*f1W)/(f1W*f2T-f2W*f1T);
      Wg = Wg + errorW;
      Tg = Tg + errorT;
      count++;
      if (count > 50){
	Wfound=-1;
	Tfound=-1;
	cout << "Newton method  doesn't converge" << endl;
	cout << "call for bisection method" << endl;
      }
    }
    Wfound=Wg;
    Tfound=Tg;
    findW();
  }

  void findW (){
    double c;
    double a=1;
    double b=sqrt(1+S*S/(RhoStar*RhoStar));
    double Ta, Tb, Tc;
    int count=0;
    Ta=findT(a);
    Tb=findT(b);
    while (fabs(a-b)>precision){
      c=(a+b)/2;
      Tc=findT(c);
      if (Equation1(a, Ta)*Equation1(c, Tc)<=0){
	b=c;
	Tb=Tc;
      }
      if (Equation1(b, Tb)*Equation1(c, Tc)<=0){
	a=c;
	Ta=Tc;
      }
      count++;
      if (count>100){
	cout <<"solution does not converge" <<endl;
	break;
      }
    }
    c=(a+b)/2;
    Wfound=c;
    Tfound=findT(Wfound);
    if (count>100){
      Wfound=-1;
      Tfound=-1;
    }
  }
  
  double findT(double W){
    double c;
    double a=-1.0;
    double b=1.;
    int count =0;
    if (!polynomial){
      if(Equation2(W, a)*Equation2(W, b)>0)
	{
	  cout << "T can't be found" << endl;
	  return -1.;
	}
      while (fabs(a-b)>precision){
	c=(a+b)/2;
	if (Equation2(W, a)*Equation2(W, c)<=0){
	  b=c;
	}
	if (Equation2(W, b)*Equation2(W, c)<=0){
	  a=c;
	}
	count++;
	if (count>100){
	  cout << "solution does not converge" << endl;
	  return -1;
	}
      }
      c=(a+b)/2;
      return c;
    }else{
      double P=-(W-1-Tau/RhoStar)/(2*W*W/(Detg*RhoStar)-Detg/RhoStar);
      double T=(P-100*(RhoStar/(Detg*W))*(RhoStar/(Detg*W)))/(RhoStar/(Detg*W));
      return T;
    }
  }
  
  vector<double> MakeInversion(double Winit, double Tinit){
    auto t1 = std::chrono::high_resolution_clock::now();
    /*if want to use bisection method only call
     findW() function instead of findWT(Winit, Tinit)*/
    //findW();
    findWT(Winit, Tinit);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << "duration is " << duration << " ms" <<endl;
    vector<double> result (2);
    result[0]=Wfound;
    result[1]=Tfound;
    return result;
  }
};


int main(){
  //the primitive variables we need to recover
  double Detg=1, Rho0=0.0001, T=0.0002;
  vector <double> u_i= {0.7,0.1};
  //computation of conservative variables
  double W=1;
  for (int i=0; i<u_i.size(); i++)
    W+=u_i[i]*u_i[i];
  W=sqrt(W);
  cout.precision(16);
  cout << "exact W = " <<W << " T= " << T << endl;
  double RhoSt = RhoStar(Detg, W, Rho0);
  double tau = Tau(Detg, W, Rho0, T);
  vector <double> S_k= Sk(Rho0, T, W, Detg, u_i);


  /* if you wanna check the bisection method for 
     the temperature evaluation - change last argument to 0*/
  Inversion A(RhoSt, tau, S_k, 1);
  vector<double> WT=A.MakeInversion(2, 0.0005);


  if ((WT[0]>=1)&&(WT[1]>=0)){
      cout << " we found W=" << WT[0] << " T=" << WT[1] << endl;
  } else{
    cout << "inversion was failed" << endl;
  }
  return 0;
}
