/* Measure angle of e1, e2, and e3 equilibrium triad arms.
   Compile: g++ -std=c++11 calc_triad_angle.cpp
   Execute: ./a.out 
*/

#include <algorithm> // for std::transform,find,reverse
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <string>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>
#include <list>
#include <map>
#include <vector>
#include <string>

using namespace std;

double rad2deg = 180/3.1415926535;
int getNframes(string ifn);
int readTriad(string ifn, double ***e_v); 
double dotprod(double a[], double b[]);
double getAngle(double *bvec1, double *bvec2, char flag);

/****************************************************************/
int main(int argc, char *argv[])
{
  int j,i,nfrm,nfr0,nfr1;
  double **phi, ***e_va, ***e_vb; 
  string if0; // input file name

  vector<int> pc; pc.push_back(2);

  if0="./data_vab/va.pdb";
  nfrm=getNframes(if0)/6; // divide by 6 b/c of 6 entries per frame

  // create variables to store triads 
  e_va=new double**[nfrm]; e_vb=new double**[nfrm]; phi=new double*[3];
  for (i=0;i<nfrm;i++) {
    e_va[i]=new double*[3];e_vb[i]=new double*[3];
    for (j=0;j<3;j++) {
      e_va[i][j]=new double[3];      e_vb[i][j]=new double[3];
    }
  } //  for (i=0;i<nfrm;i++) {
  for (i=0;i<3;i++) phi[i]=new double[nfrm]; 
  
  // Read data
  if0="./data_vab/va.pdb";       nfr0=readTriad(if0,e_va);
  if0="./data_vab/vb.pdb";       nfr1=readTriad(if0,e_vb);
  
  assert(nfrm==nfr0);  assert(nfrm==nfr1); // sanity check
  
  // Get angle per frame
  for (i=0;i<nfrm;i++) {
    for (j=0;j<3;j++) // 3 triad arms
      phi[j][i]=getAngle(e_va[i][j],e_vb[i][j],3);
  }

  ofstream of1, of2, of3; // print each triad angle to separate file
  of1.open ("phi1.txt");   of2.open ("phi2.txt");   of3.open ("phi3.txt"); 

  for (i=0;i<nfrm;i++) of1 <<" "<<phi[0][i]*rad2deg<<endl;
  for (i=0;i<nfrm;i++) of2 <<" "<<phi[1][i]*rad2deg<<endl;
  for (i=0;i<nfrm;i++) of3 <<" "<<phi[2][i]*rad2deg<<endl;
  
  of1.close(); of2.close(); of3.close(); 
  return 0;
}

/************************************************************************/
int getNframes(string ifn)
{
  int nfrm=0;
  string line;
  ifstream IFN (ifn);
  if (IFN.is_open()) {
    while (getline(IFN,line)) ++nfrm;
  }
  else cout<<"File not opened."<<endl; 

  IFN.close();
  return nfrm; 
}

/************************************************************************/
int readTriad(string ifn, double ***e_v)
{ /* Read file and return e* triads as such:
     triad[0-2]=e_va, triad[3-5]=e_vb.
     - entries: cm, e1, e2, e3
  */

  int iresi,j,nfrm=0;
  double idum, etri[3],tvec[3],mag,cm[3];
  string sdum; 
  
  string line;
  ifstream IFN (ifn);
  if (IFN.is_open()) {
    while (getline(IFN,line))
      {
	std::stringstream ss(line);
	//if (nfrm==5) break; 
	if (line=="MODEL") continue;
	else if (line=="ENDMDL") nfrm++;
	if (ss >> sdum >> iresi >> sdum >> sdum >> idum >> 
	    etri[0] >> etri[1]>> etri[2] >> idum >> idum >> sdum) {
	  if (iresi==1) { // cm coordinate
	    for (j=0;j<3;j++) cm[j]=etri[j];
	  }
	  else { // non-cm, calc unit vect b/en cm and triad
	    for (j=0;j<3;j++)  tvec[j]=etri[j]-cm[j];
	    mag=sqrt(tvec[0]*tvec[0]+tvec[1]*tvec[1]+tvec[2]*tvec[2]);
	    for (j=0;j<3;j++) 	{
	      e_v[nfrm][iresi-2][j]=tvec[j]/mag;
	    }
	  }
	}
      }
  } //   if (IFN.is_open()) {
  else cout<<"Unable to open file: "<<ifn<<endl;
  
  IFN.close();
  
  return nfrm; 
}

/************************************************************************/
double dotprod(double a[], double b[])
{ // Calculates \vec a\cdot \vec b in 3D

  double r=0.0;  
  for (int i=0;i<3;i++) r+=a[i]*b[i];
  return r;
}

/************************************************************************/
double getAngle(double *bvec1, double *bvec2, char flag)
{ /* Find angle between bvec1 & bvec2. 
     flag={0,1}: 2-Dim angle. 
        flag=0: Range is (0,PI)
        flag=1: Range is (-PI,PI). The sign is defined in the
           counterclockwise manner from bvec1 to bvec2.  
     flag=3: 3-dim angle
 */
  
  int ndim;
  double r1, r2, r3, phi;
  double TOL = 1.0e-5;
    
  ndim=(flag<3)?2:3;
  r1=dotprod(bvec1,bvec1);  r2=dotprod(bvec2,bvec2);
  r1=sqrt(r1); r2=sqrt(r2);
  r3=dotprod(bvec1,bvec2)/(r1*r2);
  r3=(r3>1.0)?(1.0-TOL):r3; // take care of round off error
  r3=(r3<-1.0)?(-1.0+TOL):r3;
  phi=acos(r3);
  if (flag==1) { // get the sign 
    r3=bvec1[0]*bvec2[1]-bvec1[1]*bvec2[0]; // cross product
    if (r3<0) phi*=-1.0;
  }
  return phi;
}
