#include "get_bond.h"

/************************************************************************/
contact::contact() { // constructor
  int i,j;
  //traj=new char*[maxbond];
  //  occ=new double[maxbond];
  for (i=0;i<maxbond;i++) {
    resi1[i]=resi2[i]=0; //traj[i]=new char[maxframe];
    //for (j=0;j<maxframe;j++) traj[i][j]=0;
  }
  nframe=nbond=0;
}

/************************************************************************/
double contact::calculate_occupancy(int i0, int n1, int n2)
{ /* calculate occupancy of bond i0 in interval [n1,n2) */
  int j;
  double rdum=0.;
  for (j=n1;j<n2;j++) rdum+=(double)traj[i0][j];
  return rdum/(double)(n2-n1);
}

/************************************************************************/
void contact::get_args(string cfname)
{
  ifstream ff(cfname.c_str());
  string sdum1; char cdum[256];
  //ifname=ofname=occ_name=btype="";
  ifname=ofname=btype="";
  btype="hbond"; // default
  npad=-1; lcut=hcut=dt=-1.;
  trajname="none"; // default value
  frm_ini=0; frm_fin=-1; stride=1; // default
  osele="";

  while (!ff.eof()) {
    ff>> sdum1;
    if (sdum1.find('#')==0) { // !=std::string::npos) {
      ff.getline(cdum,256); // comment. ignore the rest
    }
    else if (strcmp(sdum1.c_str(),"ifname")==0) ff>>ifname;
    else if (strcmp(sdum1.c_str(),"ofname")==0) ff>>ofname;
    else if (strcmp(sdum1.c_str(),"read_bond_data")==0) ff>>trajname;
    else if (strcmp(sdum1.c_str(),"btype")==0) ff>>btype;
    else if (strcmp(sdum1.c_str(),"npad")==0) ff>>npad;
    else if (strcmp(sdum1.c_str(),"lcut")==0) ff>>lcut;
    else if (strcmp(sdum1.c_str(),"hcut")==0) ff>>hcut;
    else if (strcmp(sdum1.c_str(),"segname")==0) ff>>segname;
    else if (strcmp(sdum1.c_str(),"dt")==0) ff>>dt;
    else if (strcmp(sdum1.c_str(),"frm_ini")==0) ff>>frm_ini;
    else if (strcmp(sdum1.c_str(),"frm_fin")==0) ff>>frm_fin;
    else if (strcmp(sdum1.c_str(),"stride")==0) ff>>stride;
    else if (strcmp(sdum1.c_str(),"osele")==0) ff>>osele;
    else {
      cout<< sdum1<<": Unrecognized option in " << cfname<<endl;
    }
  }
  if ((strlen(ifname.c_str())==0)&&(trajname=="none")) {
    cout<<"ERROR: No ifname in "<<cfname<<endl; exit(-1);
  }
  if  (strlen(ofname.c_str())==0) {
    cout<<"ERROR: No ofname in "<<cfname<<endl; exit(-1);
  }

  hcut0=0.5*hcut;
  return;
}
    
/************************************************************************/
void contact::get_occupancy() 
{ /* Calculate occupancy and lifetime */
  int i,j, nf0=(frm_fin-frm_ini)/stride;
  double rdum;
  for (i=0;i<nbond;i++) {
    rdum=0.;
    for (j=frm_ini;j<frm_fin;j+=stride) rdum+=(double)traj[i][j];
    rdum/=(double)nf0;
    occ_bond.insert(pair<double,int>(rdum,i));
    bond_occ.insert(pair<int,double>(i,rdum));
  }  
}

/************************************************************************/
void contact::get_occupancy(int frm0, int frm1) 
{ /* Calculate occupancy between frame frm0 and frm1. */
  int i,j, nf0=(frm1-frm0); // /stride; (stride unused)
  double rdum;
  for (i=0;i<nbond;i++) {
    rdum=0.;
    for (j=frm0;j<frm1;j+=1) rdum+=(double)traj[i][j];
    rdum/=(double)nf0;
    occ_bond.insert(pair<double,int>(rdum,i));
    bond_occ.insert(pair<int,double>(i,rdum));
  }  
}

/************************************************************************/
void contact::get_occ_traj() 
{ /* Find occupancy trajectory  */
  int i,j,npad0=npad/2;
  double occ;

  if (nframe < (2*npad)) {
    cout<< "Too small number of frames. Reduce npad in get_bond.cpp."<<endl;
    exit(-1);
  }
  cout <<"Calculating occupancy trajectory with npad= "<<npad<<endl;
  for (i=0;i<nbond;i++) {
    // first check if first and last npad frames indicate transition
    for (j=0;j<(nframe-npad);j++) {
      occ = calculate_occupancy(i,j,j+npad);
      occ_traj[i][j+npad0]=occ;
    }
  }
}

/************************************************************************/
void contact::get_transition() 
{ /* Find bond transition. Algorithm:
    1) Use occ_traj for operations below.
    2) Use occupancy for first and last npad/2 points and check for 
       transitions:
    2a) occ<lcut, occ>hcut at beginning/end respectively: find bond formation
    2b) occ<hcut, occ>lcut at beginning/end respectively: find bond breakage.
    2c) Otherwise: No transition. skip to the next bond.
    3) Find transition time for 2a) and 2b):
      In the case of bond breakage, reverse occ_traj and save into otrj.
      3a) Find first frame t0 where occ>lcut. Transition time is between
          this time t0 and end time nframe0=nframe-npad/2
      3b) Find avg occ between t0 and nframe0, called occ2
      3c) Increase t from t0, and find first frame t1 where occ>(0.5*occ2)
      3d) Find avg & std of occ between t1 and nframe0, called occ2
      3e) In the case of bond breakage, t1=nframe0-1-t1
          store t1, occ1, std(occ1).
  */
  int i,j,k,npad0=npad/2, nframe0=nframe-npad0, t0,t1,t2,t3,flag;
  int flag_t; // flag for transition type
  double occ,occ0,occ1,occ2,sig0,sig1,sig2;
  double rdum1,rdum2,rdum3,rdum4, *otrj=new double[nframe];
  cout<<"Finding bond formation/breakage events."<<endl;
  for (i=0;i<nbond;i++) {
    getavg(&occ0, &sig0, &occ_traj[i][npad0], npad0);
    getavg(&occ1, &sig1, &occ_traj[i][nframe-npad], npad0);
    // skip if no clear transition
    if ((occ0<=hcut)&&(occ0>=lcut)) continue; // intermed begin
    if ((occ1<=hcut)&&(occ1>=lcut)) continue; // intermed end
    if ((occ0<lcut)&&(occ1<lcut)) continue; // both low occ
    if ((occ0>hcut)&&(occ1>hcut)) continue; // both high occ
    // initial criterion passed. 
    if (occ0 <lcut) { // search bond formation
      assert(occ1 > hcut);
      for (j=0;j<nframe;j++) otrj[j]=occ_traj[i][j];
      flag_t=0;
    }
    else { // search bond breakage (reverse search)
      assert(occ0>hcut); assert(occ1<lcut);
      for (j=0;j<nframe;j++) otrj[nframe-j-1]=occ_traj[i][j];
      flag_t=1;
    } // else { // search bond breakage
    rdum1=(flag_t==0)?occ1:occ0; rdum1*=0.5;
    flag=0;
    for (j=npad0;j<nframe0;j++) { // 3a)
      if (otrj[j]>lcut) {t0=j; flag=1; break;} 
    }
    assert(flag==1); flag=0;
    getavg(&occ2,&sig2,&otrj[t0],(nframe0-t0)); // 3b)
    rdum2=0.5*occ2;
    for (j=t0;j<nframe0;j++) { // 3c)
      if (otrj[j]>rdum2) { t1=j; flag=1; break;}
    }
    assert(flag==1);
    getavg(&occ2,&sig2,&otrj[t1],(nframe0-t1)); // 3d)
    j=(flag_t==0)?t1:(nframe-1-t1);
    transition.insert(pair<int,int>(j,i));
    occ2=(flag_t==0)?occ2:-1.*occ2; // neg occupancy for bond breakage
    trans_occ.insert(pair<int,double>(i,occ2)); 
    trans_occ_std.insert(pair<int,double>(i,sig2)); 
  } // for (i=0;i<nbond;i++) {
}

/************************************************************************/
void contact::read_data()
{
  ifstream ff(ifname.c_str());
  string line,sdum,sdum1,sdum2,sdum3,sdum4, linedata[12];
  int i,j,idum,idum1,idum2,flag;
  size_t found2;
  char **traj_tmp=new char*[maxbond];
  for (i=0;i<maxbond;i++) {
    traj_tmp[i]=new char[maxframe];
    for (j=0;j<maxframe;j++) traj_tmp[i][j]=0;
  }
  
  while (!ff.eof()) {
    /* specify which phrase to use as a beginning of a new block in 
       charmm output */
    if (strcmp(btype.c_str(),"hbond")==0) 
      found2=line.find("I-atom");
    else if (strcmp(btype.c_str(),"nonpolar")==0) 
      found2=line.find("MIN DISTANCE");
    else { cout<<"ERROR: Wrong bond type!"<<endl; exit(-1);}

    if (found2!=std::string::npos) {
      getline(ff,line); 
      // read one more line in case of hbond
      if (strcmp(btype.c_str(),"hbond")==0) getline(ff,line); 
      std::size_t found3=line.find(segname.c_str());
      while (found3!=std::string::npos) { // read in distance block
	stringstream ss(line);
	for (i=0;i<8;i++) ss>>linedata[i];
	if (strcmp(btype.c_str(),"hbond")==0) {
	  sdum1=linedata[1]; sdum2=linedata[6]; // resname
	  sdum3=linedata[0]; sdum4=linedata[5]; // segid
	  idum1=atoi(linedata[2].c_str());  // resids for h-bonds
	  idum2=atoi(linedata[7].c_str());
	}
	else if (strcmp(btype.c_str(),"nonpolar")==0) {
	  for (i=8;i<12;i++) ss>>linedata[i];
	  sdum1=linedata[2]; sdum2=linedata[8]; // resname
	  sdum3=linedata[1]; sdum4=linedata[7]; // segid
	  idum1=atoi(linedata[3].c_str());  // resids
	  idum2=atoi(linedata[9].c_str());
	}
	if ((sdum3==sdum4)&&(idum1 == idum2)) { //skip the same resid pair
	  getline(ff,line); found3=line.find(segname.c_str());
	  continue;
	}
	else if ((sdum3==sdum4)&&(idum1>idum2)) { // order resi w/in same segi
	  idum=idum1; sdum=sdum1;
	  idum1=idum2; sdum1=sdum2; idum2=idum; sdum2=sdum;
	}
	else {}
	flag = 0; // flag for existing bond
	for (i=0;i<nbond;i++) {
	  if ((sdum3==segi1[i])&&(idum1 == resi1[i])) {
	    if ((sdum4==segi2[i])&&(idum2 == resi2[i])) {
	      // existing bond
	      if (traj_tmp[i][nframe] == 0) { // avoid double counting
		traj_tmp[i][nframe] = 1; // bond i formed at nframe
	      }
	      flag=1; break;
	    }
	  }
	} // for (i=0;i<nbond;i++) {
	if (flag == 0) { // new bond
	  segi1[nbond]=sdum3; resn1[nbond]=sdum1; resi1[nbond]=idum1;
	  segi2[nbond]=sdum4; resn2[nbond]=sdum2; resi2[nbond]=idum2;
	  traj_tmp[nbond][nframe]=1;
	  ++nbond;
	} // if (flag == 0) { // new bond
	getline(ff,line); found3=line.find(segname.c_str());
      } //while (found3!=std::string::npos) {
      ++nframe;
      if (nframe%500==0) cout << nframe<<" frames read"<<endl;
    } // if (found2!=std::string::npos) {
    getline(ff,line); 
  }
  cout << nframe<<" frames read in total."<<endl;

  traj=new char*[nbond];
  if (frm_fin==-1) frm_fin=nframe; // default
  if (frm_ini<0) frm_ini=0;
  
  for (i=0;i<nbond;i++) {
    traj[i]=new char[frm_fin-frm_ini];
    for (j=frm_ini;j<frm_fin;j++) traj[i][j-frm_ini]=traj_tmp[i][j];
  }
  nframe=frm_fin-frm_ini;

  // initialize occ_traj
  occ_traj=new double*[nbond];
  for (i=0;i<nbond;i++) occ_traj[i]=new double[nframe];

  for (i=0;i<maxbond;i++) delete[] traj_tmp[i];
  delete [] traj_tmp;
};

/************************************************************************/
void contact::read_header(ifstream& ff)
{ // read file header from ff
  string sdum; 
  ff >> sdum>>sdum>>sdum>>ifname; // original input filename
  ff >> sdum >> sdum >> dt >> sdum >> sdum >> sdum >> sdum >> sdum;
  getline(ff,sdum); 
  getline(ff,sdum);
  //ff >> sdum >> sdum >> dt >> sdum >> lcut >> sdum >> hcut;
  ff >> sdum >> sdum >> sdum >> sdum >> sdum >> nframe >> sdum;
  getline(ff,sdum); // read off prev line
  getline(ff,sdum);
  ff >> sdum >> sdum >> sdum >> sdum >> nbond;
}

/************************************************************************/
void contact::read_bond_data()
{ // read trajectory from existing file
  int i,j, idum, nf0;
  double rdum;
  ifstream ff(trajname.c_str());
  string line,sdum,sdum1,sdum2, linedata[8];

  read_header(ff);
  traj=new char*[nbond];
  if (frm_fin==-1) frm_fin=nframe;
  
  for (i=0;i<nbond;i++) traj[i]=new char[nframe];
  
  // initialize occ_traj
  cout<< "Reading bond trajectory from "<<trajname<<endl;
  occ_traj=new double*[nbond];
  for (i=0;i<nbond;i++) occ_traj[i]=new double[nframe];
  getline(ff,line); // end of previous line
  getline(ff,line); // "# Bond list: "
  for (i=0;i<nbond;i++) {
    ff >> sdum >> sdum >> segi1[i] >> resn1[i] >> resi1[i]
       >> segi2[i] >> resn2[i] >> resi2[i];
  }
  getline(ff,line); // end of previous line
  getline(ff,line); // blank line
  getline(ff,line); // "#frame ..."
  for (j=0;j<nframe;j++) { // read all frames
    getline(ff,sdum);
    stringstream ss(sdum);
    ss >> rdum; // read time
    for (i=0;i<nbond;i++) { ss >> idum; traj[i][j]=(char)idum;}
    //for (i=0;i<nbond;i++) { ss >> idum; traj[i][j-frm_ini]=(char)idum;}
  }
}

/************************************************************************/
void contact::write_data()
{
  string ofn=ofname+".dat";
  ofstream ff(ofn.c_str());
  double t;
  int i,j,idum=0;
  idum=frm_ini; // if (tmode==2) idum=nframe;

  cout<< "Writing trajectory to "<<ofn<<endl;
  write_header(ff);
  ff << "# Bond list: "<<endl;
  for (i=0;i<nbond;i++) {
    ff << "# "<<setw(4)<<i<<"  "<<setw(4)<<segi1[i]<<" "
       <<setw(3)<<resn1[i]<<" "<<setw(4)<<setfill('0')<<resi1[i] <<" "
       <<setfill(' ')
       <<setw(4)<<segi2[i]<<" "<<setw(3)<<resn2[i]<<" "
       <<setw(4)<<setfill('0')<<resi2[i]<<setfill(' ')<<endl;
  }
  ff <<endl<<"#frame ";
  for (i=0;i<nbond;i++) ff <<setw(2)<<i <<setw(1)<<" ";
  ff<<endl;
  for (i=0;i<nframe;i++) {
    t= (double)(i+idum)*(double)dt;
    ff <<fixed<<setw(9)<<setprecision(3)<<t<<" ";
    for (j=0;j<nbond;j++) {
      ff<<setw(2)<< (int)traj[j][i]<<setw(1)<<" ";
    }
    ff << endl;
  }
}

/************************************************************************/
void contact::write_header(ofstream& ff)
{
  if (trajname=="none") ff << "# input filename: "  << ifname <<endl;
  else ff << "# input bond data: "  << trajname <<endl;
  ff << "# dt(ns)= "<< setw(5)<<setprecision(3)<< dt
     << "   npad= "<< setw(6)<< npad
     << "   lcut= "<< setw(5)<<setprecision(3)<< lcut
     << "   hcut= "<< setw(5)<<setprecision(3)<< hcut
     <<endl;
  ff<<"# npad<0: no local occ traj, [lh]cut<0.: no transition calculated."
    <<endl;
  
  ff << "# number of frames read: "  << nframe <<" "<<setw(8)
     <<dt*(double)nframe << "ns"<<endl;
  ff << "# frame range for analysis:  "<<setw(6)<<frm_ini<<" "
     << setw(6) << frm_fin<<endl;
  ff << "# number of bonds: " << nbond <<endl;
}
  
/************************************************************************/
void contact::write_occupancy()
{
  string ofn3=ofname+"_occ.dat";
  ofstream ff3(ofn3.c_str());
  int b, iframe, i;
  double occ,sig,rdum;
  multimap<double,int>::reverse_iterator rt;
  int idum=0;
  idum=frm_ini;

  cout << "Writing bond occupancy to "<<ofn3<<endl;
  write_header(ff3);
  ff3 << "# ind: index, bnum: bond number" <<endl;
  ff3 << "# Bond type and occupancy (ordered w/ occupancy): "<<endl;
  for (rt=occ_bond.rbegin();rt!=occ_bond.rend();rt++) {
    i=rt->second; rdum=rt->first;
    //ff3 << "# "<<setw(4)<<i<<"  "<<setw(4)<<segi1[i]<<" "
    ff3 << setw(4)<<i<<"  "<<setw(4)<<segi1[i]<<" "
	<<setw(3)<<resn1[i]<<" "<<setw(4)<<setfill('0')<<resi1[i] <<" "
        <<setfill(' ')
	<<setw(4)<<segi2[i]<<" "<<setw(3)<<resn2[i]<<" "
	<<setw(4)<<setfill('0')<<resi2[i]<<" "
	<<setfill(' ')<<setw(8)<<setprecision(5)<<fixed<<rdum<<endl;
  }

}

/************************************************************************/
void contact::write_occ_traj()
{
  int npad0=npad/2;
  string ofn1=ofname+"_occ_traj.dat";
  ofstream ff(ofn1.c_str());
  double rdum,t;  int i,j;
  multimap<double,int>::reverse_iterator it;
  int idum=0;
  idum=frm_ini; 
  
  cout<< "Writing occupancy trajectory to "<<ofn1<<endl;
  write_header(ff);
  ff << "# Bond type: "<<endl;
  for (i=0;i<nbond;i++) {
    ff << "# "<<setw(4)<<i<<"  "<<setw(4)<<segi1[i]<<" "
       <<setw(3)<<resn1[i]<<" "<<setw(4)<<setfill('0')<<resi1[i] <<" "
       <<setfill(' ')
       <<setw(4)<<segi2[i]<<" "<<setw(3)<<resn2[i]<<" "
       <<setw(4)<<setfill('0')<<resi2[i]<<setfill(' ')<<endl;
  }
  ff <<endl<<"#frame ";
  for (i=0;i<nbond;i++) ff <<setw(2)<<i <<setw(1)<<" ";
  ff<<endl;
  //  for (i=npad0;i<nframe-npad0;i++) {
  for (i=0;i<nframe;i++) {
    t= (double)(i+idum)*(double)dt;
    ff <<fixed<<setw(9)<<setprecision(3)<<t<<" ";
    for (j=0;j<nbond;j++) {
      ff<<setw(5)<< occ_traj[j][i]<<setw(1)<<" ";
    }
    ff << endl;
  }
}

/************************************************************************/
void contact::write_transition()
{
  string ofn=ofname+".dat";
  string ofn1=ofname+"_break.dat",ofn2=ofname+"_form.dat";
  ofstream ff1(ofn1.c_str());  ofstream ff2(ofn2.c_str());
  multimap<int,int>::iterator it;
  multimap<int,double>::iterator jt,kt;
  int b, iframe, i,j;
  double occ,sig,rdum;
  int idum=0;

  cout << "Writing bond transition data to "<<ofname<<"_{break,form}.dat"<<endl;
  write_header(ff1);  write_header(ff2);
  //ff1 << "# ind: index, bnum: bond number" <<endl;
  ff1 << "# occ,std(occ): occupancy and std while the bond is formed "<<endl;
  ff1 <<endl;
  ff1 << "# Broken bonds:"<<endl;
  ff1<< "# bnum                            t(frame)    t(ns)  occ    std(occ))"<<endl;

  //ff2 << "# ind: index, bnum: bond number "<<ofn <<endl;
  ff2 << "# occ,std(occ): occupancy and std while the bond is formed "<<endl;
  ff2 <<endl;
  ff2 << "# Formed bonds:"<<endl;
  ff2<< "# bnum                            t(frame)    t(ns)  occ    std(occ))"<<endl;

  i=j=0;
  for (it=transition.begin();it!=transition.end();it++) {
    iframe=it->first; b=it->second; 
    jt=trans_occ.find(b); occ=jt->second;
    kt=trans_occ_std.find(b); sig=kt->second;
    if (fabs(occ)<hcut) continue; // skip low-occupancy bonds
    if (occ <0.) {
      //ff1 << "# "<<setw(4)<<i<<"  "<<setw(4)<<b<<" "
      ff1 << setw(4)<<b<<" "
	  <<setw(4)<<segi1[b]<<" "
	  <<setw(3)<<resn1[b]<<" "<<setw(4)<<resi1[b] <<" "
	  <<setw(4)<<segi2[b]<<" "<<setw(3)<<resn2[b]<<" "
	  <<setw(4)<<resi2[b]
	  <<" " <<setw(6)
	  <<setfill('0')<<iframe<<setfill(' ')
	  << " " <<setw(9)<<setprecision(3)<<fixed<< dt*(double)iframe
	  <<" "<<setw(8)<<setprecision(5)<<fixed<< -1.0*occ
	  <<" "<<sig<<endl;
      ++i;
    }
    else {
      //ff2 << "# "<<setw(4)<<i<<"  "<<setw(4)<<b<<" "
      ff2 <<setw(4)<<b<<" "
	  <<setw(4)<<segi1[b]<<" "
	  <<setw(3)<<resn1[b]<<" "<<setw(4)<<resi1[b] <<" "
	  <<setw(4)<<segi2[b]<<" "<<setw(3)<<resn2[b]<<" "
	  <<setw(4)<<resi2[b]
	  <<" "<<setw(6)
	  <<setfill('0')<<iframe<<setfill(' ')
	  << " " <<setw(9)<<setprecision(3)<<fixed<< dt*(double)iframe
	  <<" "<<setw(8)<<setprecision(5)<<fixed<< occ
	  <<" "<<sig<<endl;
      ++j;
    }
  }

}

/************************************************************************/
void getavg(double *avg, double *sig, double ia[], int N)
/* Calculate average and s.d. of an array ia[] of size N. */
{
  int i; double a, b, tol=-1.e-12; 
  a=b=0;
  for (i=0;i<N;i++) { a+=ia[i]; b+=(ia[i]*ia[i]);}
  a/=(double)N; b/=(double)N;
  b=b-a*a; 
  if (b<0.) {assert(b>tol); b=0.;}
  *avg=a; *sig=sqrt(b);
  return;
}

/************************************************************************/
void getminmax(double *min, int *imin, double *max, int *imax, 
     double ia[], int N)
/* Finds minimum & maximum values and their potisionos of array ia of size N */
{
  int i, min0, max0;
  double rmin,rmax;
  min0=max0=0; 
  rmin=rmax=ia[0];
  for (i=0;i<N;i++) {
    if (ia[i]<rmin) {rmin=ia[i]; min0=i;}
    if (ia[i]>rmax) {rmax=ia[i]; max0=i;}
  }
  *min=rmin; *max=rmax;
  *imin=min0; *imax=max0;
}

/************************************************************************/
void contact::write_osele()
{
  int i,j,k, n0, frm0, frm1;
  string ofn;
  double t;
  int idum=0;
  idum=frm_ini; 

  //frm0=-1;
  for (i=0;i<nframe;i++) {
//    if (frm0==-1) frm0=i;
//    else frm1=i;
    
    stringstream ss;
    ss <<osele<<"/"<<ofname<<setw(5)<<setfill('0')<<i<<".str";
    ss >> ofn;
    ofstream ff(ofn.c_str());

    n0=0;
    for (j=0;j<nbond;j++) { if (occ_traj[j][i]>=ocut) ++n0;}
    ff<<"set n_"<<ofname<< " "<< n0<<endl;
    if (n0==0) {
      ff.close(); continue;
    }
    k=0;
    for (j=0;j<nbond;j++) {
      if (occ_traj[j][i]>=ocut) {
	ff<<"defi " <<ofname << k++ <<" sele segi "<<segi2[j]
	  <<" .and. resi " <<resi2[j]<<" end"<<endl;
      }
    }
    k=0;
    ff<<"defi "<<ofname <<"_tot sele ";
    for (j=0;j<nbond;j++) {
      if (occ_traj[j][i]>=ocut) {
	ff<<ofname<< k++;
	if (k<n0) ff<<" .or. ";
	else ff<<" end"<<endl;
      }
    }
    ff.close();
  }
}
