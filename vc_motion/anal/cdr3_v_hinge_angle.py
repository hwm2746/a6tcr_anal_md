# plot cdr3 distance vs hinge angle.

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import math 

#######################################################
### MODIFIABLE SECTION

fstart=24999 # first frame for analysis [fstart,end]

pwd0='../data/'
pwd1='/meas_'
pwd2='.dat'

pwd_cdr3='../../vv_motion/data/'
pwd1_cdr3='/cdr3.dat'

# Same as in hinge_angle.py
dBin=[ -1 for i in range(2)]
dBin[0]=1.3635260000000002
dBin[1]=1.1137319999999988

ifname=[ \
         ['Tab','purple','P',7,'dotted',2.5,'purple',0,'purple'],\
         ['WT0','gold','v',6,'solid',5,'k',1.0,'gold'],\
         ['WTlow','sandybrown','d',7,'solid',3.5,'k',1.5,'white'],\
         ['WThigh','darkgoldenrod','o',5,'solid',3,'k',0.7,'darkgoldenrod'],\
]

#######################################################

#######################################################
def get_data(idx):
    vc0,vc1 = [], []
    Y0,Y_cdr3=[],[]

    # v-c angle
    ifn=pwd0+ifname[idx][0]+pwd1+ifname[idx][0]+pwd2
    
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            vc0.append(values[0])  # phi_a
            vc1.append(values[1])  # phi_b

    # cdr3 distance
    ifn=pwd_cdr3+ifname[idx][0]+pwd1_cdr3
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            Y0.append(values[0])
            
    Y_cdr3=Y0[fstart:len(Y0)] # copy over to save
    vc0=vc0[fstart:len(vc0)]
    vc1=vc1[fstart:len(vc1)]

    assert(len(Y_cdr3)==len(vc0))
    
    ##############################
    # va-ca
    
    setrange=max(vc0)-min(vc0)
    nbin=math.ceil(setrange/dBin[0])
    #print("ifn:",idx," va-ca, bins=",nbin)
    
    Xm1 = [0] * nbin
    Ym1 = [0] * nbin
    Ystd1 = [0] * nbin
    cdr31 = [ [] for i in range(nbin) ] # store cdr3 dists based on angle sorting 
    
    data = np.histogram(vc0,bins=nbin) #,density=True)
    inds=np.digitize(vc0,data[1])

    for ibin in range(nbin):
        Xm1[ibin]=(data[1][ibin]+data[1][ibin+1])/2.0 
    for ibin in range(1,nbin+1): 
        for i in range(len(Y_cdr3)):
            if (inds[i]==ibin):
                cdr31[ibin-1].append(Y_cdr3[i]) 
                #if (ibin==20): print(i+fstart) # for checking
    # end for loop

    # get cdr3 avg & std per bin
    for ibin in range(nbin): # can start from 0
        if  (len(cdr31[ibin])==0): # cp prior entry to skip plotting this entry
            Ym1[ibin]=Ym1[ibin-1]
            Ystd1[ibin]=0 
        else : 
            Ym1[ibin]=sum(cdr31[ibin])/len(cdr31[ibin])
            if (len(cdr31[ibin])==1): 
                Ystd1[ibin]=0
            else :
                Ystd1[ibin]=stat.stdev(cdr31[ibin])

    #####
    # vb-cb
    
    setrange=max(vc1)-min(vc1)
    nbin=math.ceil(setrange/dBin[1])
    #print("ifn:",idx," vb-cb, bins=",nbin)
    
    Xm2 = [0] * nbin
    Ym2 = [0] * nbin
    Ystd2 = [0] * nbin
    cdr32 = [ [] for i in range(nbin) ] # store cdr3 dists based on angle sorting 
    
    data = np.histogram(vc1,bins=nbin) 
    inds=np.digitize(vc1,data[1])

    for ibin in range(nbin):
        Xm2[ibin]=(data[1][ibin]+data[1][ibin+1])/2.0 
    for ibin in range(1,nbin+1):  
        for i in range(len(Y_cdr3)):
            if (inds[i]==ibin):
                cdr32[ibin-1].append(Y_cdr3[i]) 
                #if (ibin==20): print(i+fstart) # for checking
    # end for loop

    # get cdr3 avg & std per bin
    for ibin in range(nbin): # can start from 0
        if  (len(cdr32[ibin])==0): 
            Ym2[ibin]=Ym2[ibin-1]
            Ystd2[ibin]=0 
        else : 
            Ym2[ibin]=sum(cdr32[ibin])/len(cdr32[ibin])
            if (len(cdr32[ibin])==1):
                Ystd2[ibin]=0
            else :
                Ystd2[ibin]=stat.stdev(cdr32[ibin])
                
    
    return Xm1, Ym1, Ystd1,Xm2, Ym2, Ystd2
    

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    
    nfile=len(ifname)
    X=[ [] for i in range(2)]     # pca amplitude 
    Y=[ [] for i in range(2)]     # cdr3 distance
    Ystd=[ [] for i in range(2)]  # cdr3 dist stdev 

    # subplots(nfiles, npc modes)
    fig,ax = plt.subplots(1,2,figsize=(13,5)) 
    fig.suptitle(r'CDR3 distance v. hinge angles', fontsize=12)

    for idx in range(nfile): 
        X[0], Y[0], Ystd[0], X[1], Y[1], Ystd[1]  =  get_data(idx)

        for j in range(2): 
            ax[j].plot(X[j], Y[j],color=ifname[idx][1],
                       marker=ifname[idx][2],
                       markersize=ifname[idx][3],
                       linestyle=ifname[idx][4],linewidth=ifname[idx][5],
                       markerfacecolor=ifname[idx][8],
                       markeredgewidth=ifname[idx][7],
                       markeredgecolor=ifname[idx][6])
            
            ax[j].fill_between(np.array(X[j]),
                               np.array(Y[j])-np.array(Ystd[j]),
                               np.array(Y[j])+np.array(Ystd[j]),
                               alpha=0.2,color=ifname[idx][1])
        
    ### end loop
            
    plt.minorticks_on()
    for j in range(2):
        ax[j].tick_params(axis="both",labelsize=12,length=7,width=2) 
            
        ax[j].yaxis.set_major_locator(MultipleLocator(1))
        ax[j].yaxis.set_minor_locator(MultipleLocator(0.5))
        ax[j].set_ylim([9.6,13.6]) 

        ax[j].xaxis.set_major_locator(MultipleLocator(10))        
        ax[j].xaxis.set_minor_locator(MultipleLocator(5))

    ax[0].set_xlim([58,88]) 
    ax[1].set_xlim([64,83])
    
    ax[0].set_ylabel(r'CDR3 Dist. ($\AA$)',fontsize=12)
    ax[0].set_xlabel(r'$\angle$TCR$\alpha$ ($^\circ$)',fontsize=12)
    ax[1].set_xlabel(r'$\angle$TCR$\beta$ ($^\circ$)',fontsize=12)
    
    plt.show()    
    quit() 
