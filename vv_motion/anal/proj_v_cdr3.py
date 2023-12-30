# plot cdr3, sorted by distance 
# Compare ^ with respective pc proj at that time 

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION 

pwd0='../data/'
pwd1='/traj0'
pwd2='.dat'

pcmode=['1','2','3'] # pc mode
pwd_cdr3='/cdr3.dat'

xlim=[8,16]
ylim=[-0.4,0.4]

frm_rng=[24999,-1] # frames to include. -1=last frame

ifname=[ \
         ['Vab',r'V$\alpha\beta$'],\
         ['Tab',r'T$\alpha\beta$'],\
         ['Vab-pMHC',r'V$\alpha\beta$-pMHC'],\
         ['WT0','0'],\
         ['WTlow','Low'],\
         ['WThigh','High'],\
]
#######################################################

#######################################################
def get_data(idx,jp):
    X0=[]
    Y0=[]
    cdr3=[]

    # pc projection 
    ifn=pwd0+ifname[idx][0]+pwd1+pcmode[jp]+pwd2
    
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            X0.append(values[0])

    # cdr3 distance
    ifn=pwd0+ifname[idx][0]+pwd_cdr3
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            Y0.append(values[0])


    if (frm_rng[1]==-1): b=len(Y0)
    else: b = frm_rng[1]
    cdr3=Y0[frm_rng[0]:b] # set interval

    ## to get frames of low/high cdr3 distance
    sort_index = np.argsort(cdr3)

    # sort set by pc values
    assert(len(cdr3)==len(X0)) # sanity check 
    X1,Y1=[],[]
    X1, Y1 = zip(*sorted(zip(cdr3,X0)))

    return X1,Y1

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    
    nfile=len(ifname)
    fig,ax = plt.subplots(len(pcmode),nfile,figsize=(10,5))  
    fig.suptitle(r'V$\alpha$-V$\beta$ Triad Proj. (rad.) vs. CDR3 Dist. ($\AA$)', fontsize=12) 

    for idx in range(nfile):
        for jp in range(len(pcmode)):          
            cdr3 , pcproj = get_data(idx,jp)
            g0= ax[jp][idx].hist2d(cdr3,pcproj,bins=50,cmap=plt.cm.plasma,
                                   density=True,cmin=0.0001)

            # need to hard-set xlim and ylim w/ hist2d plot 
            ax[jp][idx].set_xlim(xlim)
            ax[jp][idx].set_ylim(ylim)
            ax[jp][idx].tick_params(axis="both",labelsize=12,length=7,width=2)
            ax[jp][idx].yaxis.set_major_locator(MultipleLocator(0.3)) 
            ####
            
            # for visual
            if (idx>0):
                ax[jp][idx].set_yticklabels([])
            if (jp < len(pcmode)-1):
                ax[jp][idx].set_xticklabels([])

            ax[jp][0].set_ylabel('PC'+str(jp+1),fontsize=12)
            ax[0][idx].title.set_text(ifname[idx][1])

    # end for loop 

    plt.show()
    quit() 
