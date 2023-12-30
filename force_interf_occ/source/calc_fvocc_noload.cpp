/* calc_fvocc_noload.cpp : Measure contact occupancy over xx-ns time window. 
   Usage:
   g++ -Wall -std=c++11 calc_fvocc_noload.cpp -o fvocc_nl
   Execute:
   ./fvocc_nl name frate ui_dname resi_excl aocut iocut
   
   This is modified from calc_fvocc.cpp for the systems without force info. 
   Variables: 
   - name is system. 
   - frate is frame rate. 
   - ui_dname is domains to include in total occupancy. 
   - resi_excl excludes any specified residues from total occupancy. 
   - aocut is overall occupancy cutoff.
   - iocut is instantaneous occupancy cutoff.
*/

#include <iostream>
#include<sstream>
#include<string> 
#include <stdio.h>
#include <cmath>
#include <cassert>
#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
#include<map>
#include<bits/stdc++.h>

using namespace std;

/* Declare functions */
void readContact(string ifname, double aocut,
		 int frate,double iocut,double dt,
		 vector<string> &resi_pair,vector<vector<double>> &occ);
void splitString(string istring, vector<double> &ovect);

/**********************************************************************/
struct triplet
{
  double dt, df, docc; 
};

int main(int argc, char *argv[])
{
  int i,j,k,k0, nfrms,ifrm,npairs,fstart; 
  string ifname; // name refers to system name
  ofstream of1; // output file  
  
  double dt=0.02 ; // occ data saving rate
  double aocut,iocut;  
  
  // store data
  vector<string> dname; // all the interfaces to check.  
  vector<double> time;   // force
  vector<string> resi_pair0, resi_pair;  // list of bond selections (segname E000.. etc) 
  vector<vector<double>> occ0,occ; // presence of bond . (bond,frame=%occ)
  vector<tuple<double,double,double>> tfo; // time, force , occupancy 
  
  /* User input variables */
  assert(argc==7) ; // make sure all entries are in 
  int frate; 
  string name, ui_dname,ofname,otag;
  string resi_excl=" "; // residue to exclude. make into list later. 
  name=argv[1]; frate=atoi(argv[2]); ui_dname=argv[3]; resi_excl=argv[4];
  aocut=atof(argv[5]); // overall occupancy
  iocut=atof(argv[6]); // occupancy over given interval
  //  otag=argv[7]; // end tag for output file name. 
  
  ofname="fvocc_"+ui_dname+".txt"; 

  /******************************/  
  /* Define domains */
  
  // w/ pmhc
  if (ui_dname=="pmhc") { // all pep and mhc
    dname.push_back("pa");    dname.push_back("pb");
    dname.push_back("ma");    dname.push_back("mb");
    goto check; 
  }
  if (ui_dname=="pep") { // peptide w/ both chains 
    dname.push_back("pa");    dname.push_back("pb");
  }
  else if (ui_dname=="pa") dname.push_back(ui_dname); // just pa
  else if (ui_dname=="pb") dname.push_back(ui_dname); // just pb
  
  if (ui_dname=="mhc")  { // mhc w/ both chains
    dname.push_back("ma");    dname.push_back("mb");
  }
  else if (ui_dname=="ma") dname.push_back(ui_dname); // just pa
  else if (ui_dname=="mb") dname.push_back(ui_dname); // just pb
  
  
  // w/ tcr intra-domain
  if (ui_dname=="intra") { // all intra-domain
    dname.push_back("vab");    dname.push_back("cab");
    dname.push_back("tcra");    dname.push_back("tcrb1");
    goto check; 
  }
  // w/ tcr intra-domain
  if (ui_dname=="intranoc") { // all intra-domain, minus ca-cb
    dname.push_back("vab");    
    dname.push_back("tcra");    dname.push_back("tcrb1");
    goto check; 
  }
  if (ui_dname=="ab") { // both a-b chains
    dname.push_back("vab");    dname.push_back("cab");
  }
  else if (ui_dname=="vab") dname.push_back(ui_dname); // just vab
  else if (ui_dname=="cab") dname.push_back(ui_dname); // just cab
  
  if (ui_dname=="vc")  { // both v-c chains
    dname.push_back("tcra");    dname.push_back("tcrb1");
  }
  else if (ui_dname=="tcra") dname.push_back(ui_dname); // just va-ca
  else if (ui_dname=="tcrb1") dname.push_back(ui_dname); // just vb-cb (no fg)
  else if (ui_dname=="fg") dname.push_back(ui_dname);  // vb-fg and cb-fg 
  
 check: 
  assert(!dname.empty()) ; // sanity check
  
  /******************************/
  
  //////////////////////////////////////////////////
  /* Read/process contact data */ 

  //cout<<"Occupancy cutoffs: avg:"<<aocut<<" inst:"<< iocut<<endl; 
  
  occ.clear(); resi_pair.clear(); 
  for (i=0;i<(int)dname.size();i++) {
    cout<<"Domain added:"<<dname[i]<<endl; 
    // h-bonds
    ifname="./in_data/"+name+"/hb_"+dname[i]+".dat"; 
    readContact(ifname,aocut,frate,iocut,dt,resi_pair0,occ0); // *0 are cleared in function
    occ.insert(occ.end(), occ0.begin(), occ0.end());
    resi_pair.insert(resi_pair.end(), resi_pair0.begin(), resi_pair0.end());
    
    // nonpolar bonds
    ifname="./in_data/"+name+"/np_"+dname[i]+".dat"; 
    readContact(ifname,aocut,frate,iocut,dt,resi_pair0,occ0); 
    occ.insert(occ.end(), occ0.begin(), occ0.end());
    resi_pair.insert(resi_pair.end(), resi_pair0.begin(), resi_pair0.end());
    //    cout<<"occ.size:"<<(int)occ.size()<<endl;     
  } //   for (i=0;i<(int)dname.size();i++) {

  /******************************/
  // set-up variables for comparison
  if (occ.empty()) nfrms=1; // dummy variable to define vectors below.
  //if nfrms<0, code will exit in next step.
  else 
    nfrms=(((int)occ[0].size()-(200/dt))*dt)/frate;  // # frames based on tows saving

  //  cout<<"nfrms:"<<nfrms<<endl;
  //  cout<<"occ[0].size()"<<(int)occ[0].size()<<endl; 

  double sum[nfrms], sum_occ[nfrms]; // store sums 
  
  if (occ.empty()) {
    cout<<"WARNING: No contacts for this interface at specified occupancy."<<endl;
    cout<<"         Code exiting with empty file. "<<endl; 
    tfo.push_back(make_tuple(0,0,0));
    goto writefile; 
  }
  
  /******************************/
  // Group contact counts into sequential intervals of frate.
  // Use same fstart as loaded systems

  for (k0=0;k0<2;k0++) { // two sets (in two files) of overlapping force data
    // hard-set start frames: see ../tows/readme.txt
    if (k0==0) fstart=200; 
    else {
      if (frate==20) fstart=210;
      else if (frate==40) fstart=220;
      else assert(1); // error
    }
    
    time.clear();

    for (i=0;i<nfrms;i++) sum_occ[i]=0; 
    ifrm=-1; npairs=(int)occ.size();
  
    for (i=0;i<(int)occ.size();i++) { // loop through contacts
      //      cout<<"resi_pair:"<<resi_pair[i]<<endl;
      if (resi_excl!=" ") {
  	if (resi_pair[i].find(resi_excl) !=string::npos) { // if found
  	  // exclude this entry. subtract from total # contacts
	  //  	  cout<<"  found exclusion"<<endl; 
  	  npairs=npairs-1;
  	  continue ; // end this check (don't add this contact) and go to next 
  	}
      }
      
      ifrm=-1;           
      for (k=0;k<nfrms;k++) sum[k]=0;
      for (j=fstart/dt;j<=(int)occ[i].size();j++) { // loop starting from fstart
	//  don't count incomplete segment	
	if ((double)j*dt+frate > (int)occ[i].size())  break; 
  	if ( (double)(remainder((j*dt-fstart),frate)) == 0 ) {
	  ++ifrm;
	  time.push_back(ifrm*frate+fstart+frate/2); // midpoint of interval
	}
  	sum[ifrm]+=occ[i][j];
      } //     for (j=0;j<(int)occ[i].size();j++) {
      
      for (k=0;k<ifrm;k++) 
	sum_occ[k]+=( sum[k] / (double)(frate/dt) );
      
    } //   for (i=0;i<(int)occ.size();i++) { // loop through contacts

    for (i=0;i<ifrm;i++) 
      tfo.push_back(make_tuple(time[i],0,sum_occ[i])); 
    
  } //   for (k0=0;k0<2;k0++) { // two sets (in two files) of overlapping force data
    
  sort(tfo.begin(),tfo.end()); //force_occ.begin(),force_occ.end());
  assert(get<0>(tfo[(int)tfo.size()-1])<(int)(occ[0].size())*dt-frate/2);
  
  // Print to file
 writefile: 
  of1.open(ofname);
  for (i=0;i<(int)tfo.size();i++) 
    of1 <<get<0>(tfo[i]) <<"  "<< get<1>(tfo[i]) <<"  "<<get<2>(tfo[i]) <<endl; 

  
  of1.close();
  return 0;   
}

/**********************************************************************/
void readContact(string ifname, double aocut,
		 int frate,double iocut,double dt,
		 vector<string> &resi_pair,vector<vector<double>> &occ)
{ /* Read data in ifname, save to occ (bond,contact prescence=0 or 1).
     Excl. occupancies that don't meet aocut. 
     
     resi_pair and occ are outputs of the function.      
  */

  int i,j;
  string sdum,sdum1,sdum2; 
  vector<double> vret, vret0;
  vector<vector<double>> occ0;
  vector<string> resi_pair0;
  //  multimap<int,int> check; 
  ifstream in(ifname.c_str());
  if (!in) { cout<<"Cannot open "<< ifname<<endl; exit(EXIT_FAILURE); }

  occ.clear(); resi_pair.clear();   
  while(!in.eof()) {
    getline(in,sdum); 
    if (sdum.find("#frame ") != string::npos) { // pushback #bonds in vector
      splitString(sdum,vret);
    }
    else if (sdum.find("# ") == string::npos && !sdum.empty() ) { // occ data
      vret.clear();       splitString(sdum,vret);
      if (occ.empty()){ // should only happen once
	for (i=1;i<(int)vret.size();i++) { 
	  vret0.clear(); vret0.push_back(vret[i]);
	  occ.push_back(vret0);  // initialize
	}
      } //       if (occ.empty()){ // should only happen once
      else { // save occ data
	for (i=1;i<(int)vret.size();i++) { //i==0 is frame#
	  occ[i-1].push_back(vret[i]);  // occ counts by bonds , else
	}
      } //       else { // save occ data 
    } //     else if (sdum.find("# ") == string::npos && !sdum.empty() ) { // occ data
    else { // save string if segi sequence found 
      if (sdum.find("E00") != string::npos || sdum.find("MHC") != string::npos
    	  || sdum.find("C00") != string::npos)   resi_pair.push_back(sdum);
    } //     else { // save strings
  }
  in.close();

  /******************************/
  //Remove entries if window average is <= iocut: DON'T DO IT HERE 

  //  cout<<" incl. all:"<<(int)occ.size()<<endl; // this is # contacts       
  // Remove entries if overall average is <= aocut
  double sum, avg;
  occ0.clear();   resi_pair0.clear();  
  for (i=0;i<(int)occ.size();i++) { // bond
    sum=0;     
    for (j=0;j<(int)occ[i].size();j++) sum+=occ[i][j]; // sum occ over all frames
    avg=sum/(int)occ[i].size(); 
    if (avg>=aocut) {
      occ0.push_back(occ[i]);
      resi_pair0.push_back(resi_pair[i]);
    }
  }
  occ.clear(); occ=occ0;
  resi_pair.clear(); resi_pair=resi_pair0;
  if (occ.empty()) cout<<"WARNING: no bonds selected"<<endl;
  cout<<" after average occ cutoff:"<<(int)occ.size()<<endl; // this is # contacts     
  
  // Remove entry if inst occupancy does not reach iocut at any frate window
  occ0.clear(); resi_pair0.clear();
  for (i=0;i<(int)occ.size();i++) { // bond
    sum=0;
    for (j=200/dt;j<(int)occ[i].size();j++) { // frames, start @ 200 ns
      if ( (double)(remainder((j*dt),frate)) == 0 ) { // start new window, but first..
	// check criteria, if met, set flag
	avg=sum/(double)(frate/dt);
	if (avg>=iocut) {
	  occ0.push_back(occ[i]);
	  resi_pair0.push_back(resi_pair[i]);
      	  break; 	 
	}
	sum=0;
      } //        if ( (double)(remainder((j*dt),frate)) == 0 ) {
      sum+=occ[i][j]; // sum occ over all frames
    } //     for (j=0;j<(int)occ[i].size();j++) { // frames
  } //   for (i=0;i<(int)occ.size();i++) { // bond
  occ.clear(); occ=occ0;
  resi_pair.clear(); resi_pair=resi_pair0;
  if (occ.empty()) cout<<"WARNING: no bonds selected"<<endl;
  cout<<" after inst. occ cutoff:"<<(int)occ.size()<<endl; // this is # contacts     
    
  return; 
}


/**********************************************************************/
void splitString(string istring, vector<double> &ovect)
{ /* Split input string, return vector of doubles from entry. */

  size_t pos = 0;
  string space_delim=" "; 
  vector<string> split{};

  ovect.clear(); 
  
  while ((pos=istring.find(space_delim)) != string::npos) {
    split.push_back(istring.substr(0,pos));
    if (!istring.substr(0,pos).empty() ) {
      if (istring.substr(0,pos)!="#frame") 
	ovect.push_back(stod(istring.substr(0,pos))); 
    }
    istring.erase(0,pos+space_delim.length());    
  }
  return; 
}
