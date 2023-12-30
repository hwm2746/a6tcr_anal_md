# angle_pep.py
# get line between va cm & vb cm 
# get plane from lsqp of peptide Ca atoms, ax+by+cz+d=0
######################################################################

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
import math
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION 

fstart=9999 # to match start frame used in force_occ.py 

pwd0='../data/'
pwd1='/pep_lsqp.txt'
pwd2='../../vv_motion/data/' 

ifname=[ \
         ['WT0','gold','v','solid',0.3,'gold','k',5],\
         ['WTlow','sandybrown','d','solid',0.3,'white','k',8.5],\
         ['WThigh','darkgoldenrod','o','solid',0.3,'darkgoldenrod','k',7],\
]

##################################################

#######################################################

def read_plane(idx):
    """
    Read plane equation from charmm. 
    Need to multiply "d" (last value) by -1.
    """

    Y0=[ [] for i in range(4) ]
    # plane from lsqp
    ifn=pwd0+ifname[idx][0]+pwd1
    
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            Y0[0].append(values[0])
            Y0[1].append(values[1])
            Y0[2].append(values[2])
            Y0[3].append(values[3])

    for i in range(len(Y0[3])):
        Y0[3][i]=-1.*Y0[3][i]

    # flip matrix dimensions        
    Y1=[ [ [] for i in range(4)] for i in range(len(Y0[0]))]
    for i in range(len(Y0[0])):
        for j in range(4):
            Y1[i][j]=Y0[j][i]

    Y1=Y1[fstart:len(Y1)]
    return Y1

##########
def read_tricm(idx):
    """
    Read va, vb triad cm positions. 
    """
    cmva=[ [] for i in range(3)]
    cmvb=[ [] for i in range(3)]

    vmod=['va.pdb','vb.pdb']
    for iv in range(2):
        ifn=pwd2+ifname[idx][0]+'/'+vmod[iv]
        with open(ifn) as f:
            lines = f.readlines()
        ln = lines[1:len(lines):6]

        if (iv==0): # va
            for i in range(len(ln)):
                values = ln[i].split()
                for j in range(3):
                    cmva[j].append(values[5+j])
        else: # vb
            for i in range(len(ln)):
                values = ln[i].split()
                for j in range(3):
                    cmvb[j].append(values[5+j])

    # flip matrix dimensions
    va=[ [ [] for i in range(3)] for i in range(len(cmva[0]))]
    vb=[ [ [] for i in range(3)] for i in range(len(cmvb[0]))]
    for i in range(len(cmva[0])): # cmva and cmvb shoudl be same size
        for j in range(3):
            va[i][j]=float(cmva[j][i])
            vb[i][j]=float(cmvb[j][i])

    assert(len(va)==len(vb)) # sanity check
    va=va[fstart:len(va)]
    vb=vb[fstart:len(vb)]    
    return va,vb

##########
def get_angle(pl,pt1,pt2):
    """ 
    Get angle between line (marked by two points to get "direction vector") and plane 
    """
    dvec= [ [] for i in range(3)]
    
    # get direction vec of line made up by two points
    for i in range(3):
        dvec[i]=pt2[i]-pt1[i]

    # now get angle
    var=abs(pl[0]*dvec[0] + pl[1]*dvec[1] + pl[2]*dvec[2])
    mag1 = np.sqrt( pl[0]*pl[0] + pl[1]*pl[1] + pl[2]*pl[2])
    mag2 = np.sqrt( dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2])
    phi = var/(mag1*mag2)
    phi=np.degrees(np.arcsin(phi))
        
    return phi

##########
def get_avg(bins,data):
    """
    Get average at x bin intervals.
    """

    avg=[]
    xpos=[]
    
    ndata = round(len(data)/bins)
    for i in range(bins):
        n0=i*ndata
        n1=i*ndata+ndata
        avg.append(np.average(data[n0:n1]))
        xpos.append( (n0+n1)/2)

    return xpos, avg

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":

    nfile=len(ifname)
    phi= [ [] for i in range(nfile)]

    xavg= [ [] for i in range(nfile)]
    phiavg= [ [] for i in range(nfile)]  

    ##########
    for idx in range(nfile):
        pln=read_plane(idx)
        va,vb=read_tricm(idx)

        #print(len(pln[0])) # should = len(va[0]) = len(vb[0])
        assert(len(pln)==len(va))
        
        for i in range(len(pln)):
            phi[idx].append(get_angle(pln[i],va[i],vb[i]))

        bins=round(len(phi[idx])/2500) # divide data into 50 ns intervals
        xavg[idx],phiavg[idx]=get_avg(bins,phi[idx])
            
    # print for check, every 50 ns
    for idx in range(nfile):
        flist=np.arange(fstart,len(phi[idx]),2500)

        #for j in flist:
        #    print("   ",j,": ",phi[idx][j])

    ########## plot ##########
    fig,ax=plt.subplots(1,figsize=(6.5,5))
    fig.suptitle('Peptide angle over time',fontsize=12)
    
    for idx in range(nfile):
        xvec=list(range(len(phi[idx])))
        xvec=[ (x+fstart)/50 for x in xvec]
        ax.plot(xvec,phi[idx],color=ifname[idx][1],linewidth=0.1,alpha=ifname[idx][4])
        
        #print(stat.mean(phi[idx]))
        #print(np.std(phi[idx]))

        # plot avg
        xavg[idx]=[ (x+fstart)/50 for x in xavg[idx]]
        ax.plot(xavg[idx],phiavg[idx],color=ifname[idx][1],
                linestyle=ifname[idx][3],linewidth=3,
                marker=ifname[idx][2],markeredgecolor=ifname[idx][6],
                markerfacecolor=ifname[idx][5],markeredgewidth=0.7,
                markersize=ifname[idx][7],
                label=ifname[idx][0])


        #print(fstart,":",phi[idx][0])    
        #print("end:",phi[idx][len(phi[idx])-1])

    plt.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(250))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))

    ax.set_xlim([175,1350]) # same as combo_fv.py 
    ax.set_ylim([37,80])
    ax.tick_params(axis="both", labelsize=12,length=7,width=2)

    ax.set_ylabel(r'$\angle$peptide ($^\circ$)',fontsize=12)
    ax.set_xlabel('Time (ns)',fontsize=12) 
    
    plt.show()
    quit()

