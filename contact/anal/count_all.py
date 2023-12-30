"""count_all.py: Count number of contacts at each frame with overall
  occupancy greater than ocut1 and instantaneous occupancy greater than
  ocut.
"""

import sys
import argparse
import struct
import numpy as np

#######################################################
### MODIFIABLE SECTION 

pwd0='../data/'

ifname=[ \
         ['WT0'],\
         ['WTlow'],\
         ['WThigh'],\
         ['Vab-pMHC'],\
]

dname=['ma', 'mb', 'pa', 'pb']

# number of frames for local averaging. Do not include first and last npad
npad=40

# for statistics in write_bond()
iframe=25000  # start @ 500 ns 

#######################################################

def get_args():
    """
    Command line options.

    suffix: suffix for occupancy file, in case occupancy is for an interval.
         ex) In data_*/[hb,np]_*0_occ.dat, suffix is 0
         To consider bond occupancy for the whole traj, omit suffix

    ifn: Interface to choose (=0: TCRab-pMHC, =1: Va-Vb )
    namedef: name definitions ("namedef_*.py" without ".py")
    ofname: output image filename without extension (pdf)

    ocut1: overall occupnacy cutoff
    ocut:  instantaneous occupancy cutoff. should be >=ocut1.

    """

    parser = argparse.ArgumentParser(description=
             'Usage: count_contact.py -ocut ocut -ofn ofname')
    parser.add_argument('-ocut', metavar='ocut', type=float)
    parser.add_argument('-ocut1', metavar='ocut1', type=float, default=0.0) 
    parser.add_argument('-ofn', metavar='ofname')

    args = parser.parse_args() # args is a dictionary
    ocut = vars(args)['ocut']
    ocut1 = vars(args)['ocut1']
    ofn = vars(args)['ofn']

    return ocut,ocut1,ofn

#######################################################
def read_bond1(idx,ocut,occ_hb,occ_np,ndomain):
    """
    Read and count number of bonds 
    """

    #for i in range(ndomain):
    # Determine number of frames from dname[0]
    # (all other dname's for molecule idx0 should have the same nframe)
    ifn0=pwd0+ifname[idx][0]+'/data_'+dname[0]+'/'
    ifn=ifn0+'hb_'+dname[0]+'_occ_traj.dat' # hbond
    tmp = np.loadtxt(ifn)
    t0=tmp[:,0]
    nframe = len(tmp)
    #nframe = endframe #len(tmp)    
    
    nhb=np.zeros((ndomain,nframe),dtype=int)
    nnp=np.zeros((ndomain,nframe),dtype=int)

    for i in range(ndomain):
        idx1 = [ sublist[0] for sublist in occ_hb[i]] # bond index
        ifn0=pwd0+ifname[idx][0]+'/data_'+dname[i]+'/'
        ifn=ifn0+'hb_'+dname[i]+'_occ_traj.dat' # hbond
        tmp = np.loadtxt(ifn)
        # if nframe != endframe:        #len(tmp):
        if nframe != len(tmp):
            print(ifn0)
            print ('ERROR: number of frames mismatch')
            quit()
        n1 = len(tmp[0,:])-1 # -1: for time column
        if (n1 == 0): continue
        if len(occ_hb[i]) == 0 : continue
        for j in range(nframe):
            n0=0
            for  k in range(n1):
                if ((tmp[j,k+1] >=ocut) and ( k in idx1 )):
                    n0 = n0+1
            nhb[i,j]=n0

    for i in range(ndomain):
        idx1 = [ sublist[0] for sublist in occ_np[i]] # bond index
        ifn0=pwd0+ifname[idx][0]+'/data_'+dname[i]+'/'
        ifn=ifn0+'np_'+dname[i]+'_occ_traj.dat' # hbond
        tmp = np.loadtxt(ifn)
        if nframe != len(tmp):
            #if nframe != endframe:  #len(tmp):            
            print ('ERROR: number of frames mismatch')
            quit()
        n1 = len(tmp[0,:])-1                    
        if (n1 == 0): continue
        if len(occ_np[i]) == 0 : continue
        for j in range(nframe):
            #tmp1=tmp[j,1:] # skip time column
            n0=0
            for  k in range(n1):
                if ((tmp[j,k+1] >=ocut) and ( k in idx1 )):
                    n0 = n0+1
            #n0=len(tmp1[tmp1>=ocut])
            nnp[i,j]=n0
            
    return nhb,nnp,t0

########################################################

def read_occ(idx,ocut1):
    """
    Read occupancy-ranked bond data. Only keep bond number and occupancy 
    for bonds with occupancy >=ocut1
    """

    ndomain=len(dname)
    occ_hb=[ [] for i in range(ndomain)]
    occ_np=[ [] for i in range(ndomain)]

    for i in range(ndomain):
        
        ifn0=pwd0+ifname[idx][0]+'/data_'+dname[i]+'/'
        ifn=ifn0+'hb_'+dname[i]+'_occ.dat' # hbond

        with open(ifn, 'r') as ff:
            for line in ff:
                if line.startswith('#'): continue # skip comment
                ldata = line.split()
                if (float(ldata[-1])< ocut1): break
                occ_hb[i].append([int(ldata[0]),float(ldata[-1])])
        ff.close()
        
        ifn=ifn0+'np_'+dname[i]+'_occ.dat' # nonpolar
        with open(ifn, 'r') as ff:  
            for line in ff:
                if line.startswith('#'): continue # skip comment
                ldata = line.split()                
                if (float(ldata[-1])< ocut1): break
                occ_np[i].append([int(ldata[0]),float(ldata[-1])])
        ff.close()
    return occ_hb,occ_np,ndomain


#######################################################
def bond_stat(nhb,nnp,ndomain,t0): 
    nframe=len(t0)
    # nframe=endframe #len(t0)    
    nframe0=nframe-iframe-npad
    bstat=np.zeros((ndomain,2))

    for k in range(ndomain):
        navg=0.; nstd=0.; nframe0=0
        for i in range(iframe,nframe-npad):
            ntot=float(nhb[k,i]+nnp[k,i])
            navg=navg+ntot
            nstd=nstd+ntot*ntot
            nframe0=nframe0+1
        navg=navg/nframe0
        nstd=nstd/nframe0
        nstd=nstd-navg*navg
        nstd=np.sqrt(nstd)
        bstat[k,0]=navg; bstat[k,1]=nstd

    return bstat

#######################################################
def write_bstat(nfile,nfrm,bstat,ocut,ocut1,ofn):

    ndomain=len(dname)
    for i in range(ndomain):
        ofname=ofn+dname[i]+'.dat'
        fout=open(ofname,'w')
        fout.write('# iframe= {:5d} npad= {:3d} ocut= {:.1f}, ocut1= {:.1f}\n'.format(iframe,npad,ocut,ocut1))
        fout.write('# mol      avg     std    nfrm\n')
        for k in range(nfile):
            fout.write('{:10s} {:6.3f} {:5.3f} {:6d}\n'.format(ifname[k][0],bstat[k][i,0],bstat[k][i,1],nfrm[k]))


#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    ocut,ocut1,ofn = get_args()
    nfile=len(ifname)
    bstat=[ [] for i in range(nfile)]
    nfrm=np.zeros(nfile,dtype=int)
    for idx in range(nfile):
        occ_hb,occ_np,ndomain = read_occ(idx,ocut1)
        nhb,nnp,t0 = read_bond1(idx,ocut,occ_hb,occ_np,ndomain)
        bstat[idx]= bond_stat(nhb,nnp,ndomain,t0)
        nfrm[idx]=len(t0)

    write_bstat(nfile,nfrm,bstat,ocut,ocut1,ofn)

    quit()
