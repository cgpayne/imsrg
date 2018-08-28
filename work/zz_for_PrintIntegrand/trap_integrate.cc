//// trap_integrate.cc ////
//// takes in a function from file (x,y) and integrates via trapezoidal rule
//// by: Charlie Payne


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <sstream>


using namespace std;


int main()
{
  int N = 40001; // length of file(s), keep consistent!
  
  ifstream gtin("M0nu_integrand_GT.txt");
  if(!gtin)
  {
    cout<<"Error: with opening gtin"<<endl;
    exit(1);
  }
  ifstream fin("M0nu_integrand_F.txt");
  if(!fin)
  {
    cout<<"Error: with opening fin"<<endl;
    exit(1);
  }
  ifstream tin("M0nu_integrand_T.txt");
  if(!tin)
  {
    cout<<"Error: with opening tin"<<endl;
    exit(1);
  }
  
  double gtval[N][2];
  double fval[N][2];
  double tval[N][2];
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<2; j++)
    {
      gtin>>gtval[i][j];
      fin>>fval[i][j];
      tin>>tval[i][j];
    }
  }
  
  ofstream fout("trap_out.txt");
  if(!fout)
  {
    cout<<"Error: with opening fout"<<endl;
    exit(1);
  }
  
  double gtsum = 0;
  double fsum = 0;
  double tsum = 0;
  for (int i=0; i<N-1; i++)
  {
    gtsum += 0.5*(gtval[i+1][1] + gtval[i][1])*(gtval[i+1][0] - gtval[i][0]);
    fsum += 0.5*(fval[i+1][1] + fval[i][1])*(fval[i+1][0] - fval[i][0]);
    tsum += 0.5*(tval[i+1][1] + tval[i][1])*(tval[i+1][0] - tval[i][0]);
  }
  
  cout<<"gtsum = "<<setprecision(12)<<gtsum<<endl;
  cout<<"fsum  = "<<setprecision(12)<<fsum<<endl;
  cout<<"tsum  = "<<setprecision(12)<<tsum<<endl;
  fout<<"gtsum              = "<<setprecision(12)<<gtsum<<endl;
  fout<<"fsum               = "<<setprecision(12)<<fsum<<endl;
  fout<<"tsum               = "<<setprecision(12)<<tsum<<endl;
  
  gtin.close();
  fin.close();
  tin.close();
  fout.close();
  
  return 0;
}


//// FIN ////
