# plot cdr3 distance vs triad arm angle 

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
pwd1='/phi'
pwd2='.txt'
triarm=['1','2','3'] # which triad arm to use

pwd_cdr3='/cdr3.dat'

# Same as in va-vb_angle.py. These values are determined in va-vb_angle.py.
dBin=[ -1 for i in range(3)]
dBin[0]=1.8091999999999984
dBin[1]=1.6252899999999997
dBin[2]=1.7116100000000003

ifname=[ \
         ['Vab','blue','P',7,'solid',1.5,'blue',0,'blue'],\
         ['Tab','purple','P',7,'dotted',1.5,'purple',0,'purple'],\
         ['Vab-pMHC','red','o',6,'dashed',2,'k',0.3,'red'],\
         ['WT0','gold','v',5,'solid',4.5,'k',1.0,'gold'],\
         ['WTlow','sandybrown','d',8.5,'solid',4,'k',1.3,'white'],\
         ['WThigh','darkgoldenrod','o',5,'solid',5,'k',0.5,'darkgoldenrod'],\
]
#######################################################

#######################################################
def get_data(idx,jp,dBin):
    X0=[]
    Y0=[]
    Y_cdr3=[]

    # read triad arm agnles
    ifn=pwd0+ifname[idx][0]+'/'+pwd1+triarm[jp]+pwd2
    
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            X0.append(values[0])

    # read cdr3 distance
    ifn=pwd0+ifname[idx][0]+pwd_cdr3
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            Y0.append(values[0])

            
    Y_cdr3=Y0[fstart:len(Y0)] # copy over to save
    X0=X0[fstart:len(X0)]

    assert(len(X0)==len(Y_cdr3))

    ##############################
    setrange=max(X0)-min(X0)
    nbin=math.ceil(setrange/dBin)
    #print("ifn:",idx," arm:",jp," bins=",nbin)
    
    ##############################
    Xm = [0] * nbin
    Ym = [0] * nbin
    Ystd = [0] * nbin
    cdr3 = [ [] for i in range(nbin) ] # store cdr3 dists based on angle sorting 
    
    data = np.histogram(X0,bins=nbin) #,density=True)
    inds=np.digitize(X0,data[1])

    ##############################
    # xpos: midpoint of bin (like in vab_angle.py)

    for ibin in range(nbin):
        Xm[ibin]=(data[1][ibin]+data[1][ibin+1])/2.0 #0-1, 1-2, .. 19-20. save to 0, 1,.. 19.
    for ibin in range(1,nbin+1):  # used to check bin label
        for i in range(len(Y_cdr3)):
            if (inds[i]==ibin):
                cdr3[ibin-1].append(Y_cdr3[i]) # store to calc avg & std outside of loop
                #if (ibin==20): print(i+fstart) # for checking
    # end for loop
    

    # get cdr3 avg & std per bin
    for ibin in range(nbin): # can start from 0
        if  (len(cdr3[ibin])==0): # cp prior entry to skip plotting this entry
            Ym[ibin]=Ym[ibin-1]
            Ystd[ibin]=0 
        else : 
            Ym[ibin]=sum(cdr3[ibin])/len(cdr3[ibin])
            if (len(cdr3[ibin])==1):
                Ystd[ibin]=0
            else :
                Ystd[ibin]=stat.stdev(cdr3[ibin])

    return Xm, Ym, Ystd, cdr3
    

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    
    nfile=len(ifname)
    X=[ [] for i in range(nfile)]     # pca amplitude 
    Y=[ [] for i in range(nfile)]     # cdr3 distance
    Ystd=[ [] for i in range(nfile)]  # cdr3 dist stdev 

    # subplots(nfiles, npc modes)
    fig,ax = plt.subplots(1,3,figsize=(15,5))
    fig.suptitle(r'CDR3 distance v. triad angles', fontsize=12)
    
    for idx in range(nfile):
        for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
            cdr3=[]
            X[idx], Y[idx], Ystd[idx],cdr3 =  get_data(idx,jp,dBin[jp])

            msize=[]
            for imrk in range(len(cdr3)):
                msize.append(math.sqrt(len(cdr3[imrk])))

            ax[jp].plot(X[idx], Y[idx],
                        color=ifname[idx][1],
                        linestyle=ifname[idx][4],linewidth=ifname[idx][5],
                        marker=ifname[idx][2],markerfacecolor=ifname[idx][8],
                        markersize=ifname[idx][3], markeredgewidth=ifname[idx][7],
                        markeredgecolor=ifname[idx][6],label=ifname[idx][0])
                        
            ax[jp].fill_between(np.array(X[idx]),
                                np.array(Y[idx])-np.array(Ystd[idx]),
                                np.array(Y[idx])+np.array(Ystd[idx]),
                                alpha=0.2,color=ifname[idx][1])

            ###########
            ## END OF LOOP
    #end     for idx in range(nfile): 

    ## LABEL AXES THAT NEED TO BE SET TO BE SAME AS VAB HISTOGRAM RANGE            
    for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
        ax[jp].yaxis.set_minor_locator(MultipleLocator(.5))
        ax[jp].yaxis.set_major_locator(MultipleLocator(1.0))            

        ax[jp].set_ylim([9.6,15]) 
        ax[jp].tick_params(axis="both",labelsize=12,length=7,width=2) 
        
        if (jp>0):
            ax[jp].set_yticklabels([])

        ax[jp].xaxis.set_minor_locator(MultipleLocator(7.5))
        ax[jp].xaxis.set_major_locator(MultipleLocator(15))

        ax[jp].set_xlabel(r'$\angle$e'+str(jp+1)+' ($^\circ$)',fontsize=12)
        
    ax[0].set_ylabel(r'CDR3 Dist. ($\AA$)',fontsize=12)

    ax[0].set_xlim([143,181])
    ax[1].set_xlim([76,111])
    ax[2].set_xlim([70,103])
    
    plt.show()    
    quit() 
