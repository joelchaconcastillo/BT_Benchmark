/*
  Benchmark: "Biased Multiobjective Optimization of Decomposition Algorihm" by Hui Li et. al
  2016
*/
#include<bits/stdc++.h>
using namespace std;
long double D1(long double g, long double theta)
{
   return (g*g) + (1.0-exp(-(g*g)/theta))/5.0;
}
long double D2(long double g, long double theta)
{
   return (g*g) + (pow(fabs(g), theta)/5.0);
}
long double S1(long double x, long double gamma)
{
  return pow(fabs(x), gamma);
}
long double S2(long double x,long double gamma)
{
   if( x>=0.0 && x<0.25) return (1.0 - pow(1.0-4.0*x, gamma))/4.0;
   else if( x>=0.25 && x<0.5) return (1.0 + pow(4.0*x-1.0, gamma))/4.0;
   else if( x>=0.5 && x<0.75) return (3.0 - pow(3.0-4.0*x, gamma))/4.0;
   else if( x>=0.75 && x<=1.0) return (3.0 + pow(4.0*x-3.0, gamma))/4.0;
}
double bt1(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-10;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if((j%2)!=0) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
double bt2(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0/5.0;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D2(yj, theta);
  if(j%2) sum2 += D2(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
double bt3(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-8;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = S1(x[0], 0.02) + sum1;
   f[1] = 1.0 - sqrt(S1(x[0], 0.02)) + sum2;
}
double bt4(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-8;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = S2(x[0], 0.06) + sum1;
   f[1] = 1.0 - sqrt(S2(x[0], 0.06)) + sum2;
}
double bt5(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-10;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = (1.0 - x[0])*(1.0-x[0]*sin(8.5*M_PI*x[0])) + sum2;
}
double bt6(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-4;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - pow(x[0],0.5+((1.5*(j-1.0))/(n-1.0)));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
double bt7(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-3;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin(6.0*M_PI*x[0]);
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
   } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
long double Q(long double z)
{
  return 4.0*z*z - cos(8.0*M_PI*z)+1.0;
}
double bt8(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-3;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - pow(x[0],0.5+((1.5*(j-1.0))/(n-1.0)));
      if(!(j%2)) sum1 += Q(D1(yj, theta));
      if(j%2) sum2 += Q(D1(yj, theta));
   } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
double bt9(vector<double> &x, vector<double> &f)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, sum3=0.0, theta = 1.0e-9;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if((j%3)==0) sum1 += D1(yj, theta);
      if((j%3)==1) sum2 += D1(yj, theta);
      if((j%3)==2) sum3 += D1(yj, theta);
   } 
   f[0] = cos(0.5*x[0]*M_PI)*cos(0.5*x[1]*M_PI) + 10.0*sum1;
   f[1] = cos(0.5*x[0]*M_PI)*sin(0.5*x[1]*M_PI) + 10.0*sum2;
   f[2] = sin(0.5*x[0]*M_PI) + 10.0*sum3;
}
int main()
{
 srand(time(0));
  vector<double>f(2), x(3); 
  for(int i = 0 ; i <100000; i++)
  {
    for(auto &a:x) a = rand()/(double)RAND_MAX;
    bt3(x, f);
    for(auto &m:f)cout <<m<<" ";
    cout<<endl;
  }
  return 0;
}
