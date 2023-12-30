/* pca_dotp.cpp: Calc dot product of pca to reveal differences
   Execute:  g++ -std=c++11 pca_dotp.cpp -o dotp
             ./dotp > out_dotp.txt 
 */

#include <algorithm> 
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
void readDir(string ifn,double **vec); 
double dotprod(double a[], double b[],int N);
void setSys(vector<string> &mols); 

/************************************************************************/

// MODIFIABLE FUNCTION

void setSys(vector<string> &mols)
{
  mols.clear();

  mols.push_back("Tab");    
  mols.push_back("WT0");      
  mols.push_back("WTlow");
  mols.push_back("WThigh");

  return ; 
}

/****************************************************************/
int main(int argc, char *argv[])
{
  int i,j;
  int imol,jmol,mol_size;
  double prod,**uvec, **vvec; // vec[18 pc modes][18-dim]
  string s0,s1;
  vector<string> mol_list; // list of structures to compare
  
  uvec=new double*[18];  vvec=new double*[18];
  for (i=0;i<18;i++) {uvec[i]=new double[18]; vvec[i]=new double[18]; }

  setSys(mol_list) ; // assign mols for comparison in this function 
  mol_size=(int)mol_list.size();

  /******************************/
  /* Comparing pc between systems */
  
  for (imol=0;imol<mol_size-1;imol++) {
    for (jmol=imol+1;jmol<mol_size;jmol++) {
      s0=mol_list[imol]+"_pca.txt";   // Vab triads
      s1=mol_list[jmol]+"_pca.txt";
      
      readDir(s0,uvec);  
      readDir(s1,vvec);
      
      for (i=0;i<18;i++) {  // 18 pc modes 
      	prod=dotprod(uvec[i],vvec[i],18);  // 18-dim normalized vector
      	cout<<imol << " "<<jmol<<" "<<i<<" "<<prod<<endl;
      } 
    } //     for (jmol=imol+1;jmol<mol_size;jmol++) {
  } //   for (imol=0;imol<mol_size-1;imol++) {
  
  return 0;
}

/************************************************************************/
double dotprod(double a[], double b[], int N)
{ // Calculates \vec a\cdot \vec b in N-dimension

  double r=0.0;  
  for (int i=0;i<N;i++) r+=(a[i]*b[i]);
  return r;
}

/************************************************************************/
void readDir(string ifn,double **vec)
{
  int i,j,k,ipc,ie=0;  // ipc: pc dirs, ie: triad arms
  double mag,e1[2][3], e2[2][3], e3[2][3]; 
  string line;
  ifstream IFN (ifn);
  ipc=-1;
  if (IFN.is_open()) {
    while (getline(IFN,line))
      {
	std::stringstream ss(line);
	if (line.find("pc")!=string::npos) {ipc++; continue;}
	else if (line.find("#")!=string::npos) {continue;}
	else {
	  if (ss >> e1[ie][0] >> e2[ie][0] >> e3[ie][0] >>
	      e1[ie][1] >> e2[ie][1] >> e3[ie][1] >>
	      e1[ie][2] >> e2[ie][2] >> e3[ie][2]) {
	  }

	  for (j=0;j<3;j++) {
	    vec[ipc][ie*9+j]=e1[ie][j];
	    vec[ipc][ie*9+j+3]=e2[ie][j];
	    vec[ipc][ie*9+j+6]=e3[ie][j];
	  }
	  ie++; if (ie==2) ie=0; 	// reset
	} // 	  if (ipc>=0) { // save to vec. 
      }
  }
  else cout<<" FILE NOT FOUND."<<endl;
  assert(ipc+1==18);
  
  // unit vector
  for (i=0;i<18;i++) { // pc mode. don't need to do here b/c pc dirs
    // are already unit vecs 
    mag=0;          
    for (j=0;j<18;j++)       mag+=(vec[i][j]*vec[i][j]);  // pc dimension
    mag=sqrt(mag);
    for (j=0;j<18;j++) vec[i][j]=vec[i][j]/mag; 
  }
  
  return;
}

