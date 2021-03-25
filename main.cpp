#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

const int n=50,l=10;
const int ntraj=1000000;
double pot_coeff=5;
double harm=pot_coeff;
double histo_coeff=0.005;

double *pot_table;
double *histo;
int x=n/2;

//return the meta potential
double meta_pot(int i)
{
  double h=0;
  int d=abs(i-n/2);
  if(d>=l) h=(d-l)*(d-l)*harm;
  
  //flat outside
  if(i<n/2-l) i=n/2-l;
  if(i>n/2+l) i=n/2+l;  
  
  return histo[i]+h;
}

int main()
{
  histo=new double[n];
  pot_table=new double[n];
  memset(histo,0,sizeof(double)*n);
  
  //prepare the potential
  memset(pot_table,0,sizeof(double)*n);
  for(int i=0;i<n;i++)
    {
      int d=i-n/2;
      //pot_table[i]=pot_coeff*pow(d,2);
      pot_table[i]=pot_coeff*cos(M_PI*d/3.0);
    }
  
  int itraj=0;
  do
    {
      //extract a test variable
      int j=(x+(random()&0x2)-1+n)%n;
      
      //compute potential
      double vx=pot_table[x]+meta_pot(x);
      double vj=pot_table[j]+meta_pot(j);
      
      //compute ptrans and perform metrotest
      double ptr=exp(-(vj-vx));
      double extr=rand()/(double)RAND_MAX;
      if(extr<ptr) x=j;
      
      //add meta_pot
      histo[x]+=histo_coeff;
      
      cout<<x<<endl;
      itraj++;
    }
  while(itraj<ntraj);
  
  //write the metapotential
  ofstream meta_pot_f("meta_pot");
  ofstream pot_f("pot");
  double meta_pot_c=meta_pot(n/2);
  for(int i=0;i<n;i++)
    {
      meta_pot_f<<i<<" "<<meta_pot_c-meta_pot(i)<<endl;
      pot_f<<i<<" "<<pot_table[i]-pot_table[n/2]<<endl;
    }
  
  delete[] pot_table;
  delete[] histo;
  
  return 0;
}
