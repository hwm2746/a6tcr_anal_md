/* hamming_dist.cpp: Calculates Hamming distances of contact patterns

  Usage:
  g++ hamming_dist.cpp ftn_bond.cpp -O3 -o te1
  ./te1 input_hamming.dat
  input_hamming.dat: input configuration file, containing:

  ndata: number of contact data to be read. Must be followed by ncontact
  lines of trajname (output of get_bond.cpp)
  ofname: output filename
  frm0/frm1: first/last frame number to use for analysis
  Default: frm0=0, frm1=nfrm
*/
#define SET_EXT // Do not put extern for global vars

#include "get_bond.h"

// global vars
int ndata,nfrm,nbond_tot,nbond_ref, *hd, *hd1, frm0,frm1;
int occ_frm; // Number of initial frames to calculate occupancy
double ocut1;
contact *c;
string ifname,ofname;

void hmg_dist();
void hmg_dist1(char **bref);
int last_frame(int frm1);
void prep_occupancy(int f0, int f1);
void prep_bond_occ();
void prep_bond_ref(char **bref);
void write_hd();

/************************************************************************/
int main(int argc, char *argv[])
{
  if (argc!=2) {
    cout<<"Usage: te1 input_hamming.dat"<<endl; return -1;
  }
  int i,j,k;
  string sdum1; char cdum[256];
  char **bond_ref;
  ifname=argv[1]; ofname="";
  
  // Read input data
  ifstream ff(ifname.c_str());
  nbond_tot=0; frm0=frm1=0; ocut1=0.; occ_frm=1;
  while (!ff.eof()) {
    ff>> sdum1;
    if (sdum1.find('#')==0) { // !=std::string::npos) {
      ff.getline(cdum,256); // comment. ignore the rest
    }
    else if (strcmp(sdum1.c_str(),"ndata")==0) {
      ff>>ndata;
      c=new contact[ndata];
      for (i=0;i<ndata;i++) {
	ff>>c[i].trajname;
	c[i].read_bond_data();
	nbond_tot+=c[i].nbond;
      }
    }
    else if (sdum1=="") continue;
    else if (sdum1=="ocut") ff>>ocut1;
    else if (sdum1=="occ_frm") ff>>occ_frm;
    else if (sdum1=="ofname") ff>>ofname;
    else if (sdum1=="frm0") ff>>frm0;
    else if (sdum1=="frm1") ff>>frm1;
    else {
      cout<< sdum1<<": Unrecognized option: " << sdum1<<endl;
    }
    sdum1.clear();
  }
  if (ofname=="") {
    cout<<" ERROR: ofname missing."<<endl;
    exit(-1);
  }

  nfrm=last_frame(0);
  cout<<"nfrm= "<<nfrm<<endl;
  if (frm1==0) frm1=nfrm; // default
  hd=new int[nfrm];
  hmg_dist();

  // get occupancy for first 2500 frames (50ns)
  hd1=new int[nfrm];
  bond_ref=new char*[ndata]; // reference bond
  for (i=0;i<ndata;i++) {
    k=c[i].nbond;
    bond_ref[i]=new char[k];
    for (j=0;j<k;j++) bond_ref[i][j]=0;
  }
  
  prep_occupancy(frm0,frm0+occ_frm);
  prep_bond_ref(bond_ref); // select bonds w/ occ>0.8
  hmg_dist1(bond_ref);

  write_hd();
  
  return 0;
}

/************************************************************************/
void hmg_dist()
{ /* get Hamming distance from ifrm, to frm1 */
  int i,j,ifrm;

  for (ifrm=frm0;ifrm<frm1;++ifrm) {
    hd[ifrm]=0;
    for (i=0;i<ndata;i++) {
      for (j=0;j<c[i].nbond;j++) {
	if (c[i].traj[j][ifrm]!=c[i].traj[j][frm0]) ++hd[ifrm];
      }
    }
  }
  return;
}

/************************************************************************/
void hmg_dist1(char **bref)
{ /* get Hamming distance from bref. 

   : Increase distance ONLY when bref for a bond is 1 and the bond
   is not formed */
  int i,j,ifrm;

  for (ifrm=frm0;ifrm<frm1;++ifrm) {
    hd1[ifrm]=0;
    for (i=0;i<ndata;i++) {
      for (j=0;j<c[i].nbond;j++) {
	//if (bref[i][j]!=c[i].traj[j][ifrm]) ++hd1[ifrm];
	if ((bref[i][j]==1)&&(c[i].traj[j][ifrm]==0)) ++hd1[ifrm];
      }
    }
  }
  return;
}

/************************************************************************/
int last_frame(int frm1)
{ /* set the last frame number to work with. If frm1>0, return frm1 */
  int i, frm1_min=1e9;
  
  if (frm1>0) return frm1;
  else {
    for (i=0;i<ndata;i++) {
      if (c[i].nframe <frm1_min) frm1_min=c[i].nframe;
    }
    return frm1_min;
  }
}

/************************************************************************/
void prep_bond_ref(char **bref)
{ /* set reference bond with occupancy higher than ocut */
  int i,j;
  multimap<int,double>::iterator it;
  nbond_ref=0;
  for (i=0;i<ndata;i++) {
    for (it=c[i].bond_occ.begin();it!=c[i].bond_occ.end();++it) {
      if ((it->second)>ocut1) {
	bref[i][it->first]=1;
	++nbond_ref;
      }
      else bref[i][it->first]==0;
    }
  }
  return;
}

/************************************************************************/
void prep_occupancy(int f0, int f1)
{
  int i;
  for (i=0;i<ndata;i++) c[i].get_occupancy(f0,f1);
}
  
/************************************************************************/
void write_hd()
{ /* write Hamming distance data */
  int i,ifrm, idum=0;
  double t;
  ofstream ff(ofname.c_str());

  cout<<"Output written to "<<ofname<<endl;
  ff<< "# ndata= "<<setw(2)<<ndata<<"  nbond_tot= "<<setw(3)<<nbond_tot
    <<"  nfrm= "<<setw(6)<<nfrm<<endl;
  ff<<"# input file: "<<ifname<<endl;
  ff<<"# col0: time(ns), col1: Hammng dist from t=0"<<endl;
  ff<<"# col2: Hamming dist from contacts with occ higher than "<<ocut1
    <<" during first "<< occ_frm<<" frames"<<endl;
  ff<<"# Number of high-occ bonds for col2: "<<nbond_ref<<endl;
  
  for (ifrm=frm0;ifrm<frm1;++ifrm) {
    t= (double)(ifrm+idum)*(double)c[0].dt;
    // ff<<setw(6)<<ifrm<<" "
    ff<<fixed<<setw(9)<<setprecision(3)<<t<<" "<<setw(6)<<hd[ifrm]<<" "
      <<setw(6)<<hd1[ifrm]<<endl;
  }
  return;
}
  
