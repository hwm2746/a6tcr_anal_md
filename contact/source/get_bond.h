/* get_bond.h */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <map>
using namespace std;

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif


#ifndef CONTACT
#define CONTACT

#define maxbond 5000
#define maxframe 100000
#define ocut 0.8 // occupancy cutoff 

class contact{
public:
  int nframe,nbond,ntrans;
  int resi1[maxbond],resi2[maxbond];
  int npad,frm_ini,frm_fin,stride; // npad: number of initial/final frames for padding
  double lcut,hcut,hcut0,dt;
  string resn1[maxbond],resn2[maxbond],segi1[maxbond],segi2[maxbond];
  string ifname,ofname,btype; //,occ_name;
  string trajname,segname;
  string osele; // for analysis
  char **traj;
  double **occ_traj; // ocupancy trajectory
  multimap<double,int> occ_bond; // int: bond number  
  multimap<int,double> bond_occ; // 
  // Transition related:
  // transition: bond number and transition frame number 
  // trans_occ: bond number and occupancy before/after tran
  // trans_occ_std: std of the local avg occupancy before/after transition
  // occupancy <0: bond breakage. >0: bond formation
  multimap<int,int> transition;
  multimap<int,double> trans_occ,trans_occ_std; 
  // Member function
  contact (); // constructor
  double calculate_occupancy(int i0, int n1, int n2);
  void get_args(string cfname);
  void get_occupancy();
  void get_occupancy(int frm0, int frm1);
  void get_occ_traj();
  void get_transition();
  void read_bond_data();
  void read_data();
  void read_header(ifstream& ff);
  void write_data();
  void write_header(ofstream& ff);
  void write_occupancy();
  void write_occ_traj();
  void write_osele();
  void write_transition();  
};

/************************************************************************/
void getavg(double *avg, double *sig, double ia[], int N);
void getminmax(double *min, int *imin, double *max, int *imax, 
	       double ia[], int N);


#endif // #ifndef CONTACT
