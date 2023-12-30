# Plot differences in PC amplitude from BOC PC 

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
import math
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION

pwd0='../data/'

pwd2='/dir_cm.dat'
pwd3='/std_cm.dat'

npc=3 # how many PC modes to plot 

ifname=[ \
         ['Tab',r'T$\alpha\beta$','purple','purple',''],\
         ['WT0','0','gold','dimgrey',''],\
         ['WTlow','L','sandybrown','black',''],\
         ['WThigh','H','darkgoldenrod','black','....'],\
]

#######################################################

######################################################################
def get_std(idx):
    Y0=[]    
    ifn=pwd0+ifname[idx][0]+pwd3
    
    for line in open(ifn,'r'):
        if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            if (values[0]>npc): 
                break; 
            Y0.append(values[1])
            
    return Y0

######################################################################
def calc_mag(idx):
    Y0=[ [] for i in range(npc)] # raw data
    Ymag=[ [] for i in range(npc)]   # magnitude
    
    ifn=pwd0+ifname[idx][0]+pwd2
    pcmode=-1
    
    for line in open(ifn,'r'):
        if line.startswith("#"):
            pcmode=pcmode+1
            #values = [string(s) for s in line.split()]
            #pcmod=int(values[2])
            if (pcmode>npc): break
        else:
            values = [float(s) for s in line.split()]
            Y0[pcmode-1].append(values)

    for i in range(npc):
        for j in range(6):
            Ymag[i].append(math.sqrt(Y0[i][j][0]*Y0[i][j][0] + Y0[i][j][1]*Y0[i][j][1] + Y0[i][j][2]*Y0[i][j][2]))

    return Ymag


#######################################################
# Main routine
#######################################################
if __name__ == "__main__":

    nfile=len(ifname)
    stdcm = [ [] for i in range(nfile)]
    vecmag = [ [ [] for i in range(npc) ] for i in range(nfile)]
    pccomp = [ [ [] for i in range(6) ] for i in range(npc)]
    rat = [ [ [] for i in range(2) ] for i in range(npc) ]

    fig,ax=plt.subplots(2,3,figsize=(7,4.5))
    fig.suptitle(r'PC amplitude per region (solid=V$\alpha$ or H$\alpha$, hashmarks=V$\beta$ or H$\beta$)') 
    figsub,axsub=plt.subplots(2,3,figsize=(7,4.5)) # ha-hb and va-vb rows, pc1--3 cols
    figsub.suptitle('PC amplitude difference')
    
    # get values 
    for idx in range(nfile):
        stdcm[idx] = get_std(idx)
        vecmag[idx] = calc_mag(idx)
        
    # scale by std multiplication
    for idx in range(nfile):
        for i in range(3): # b/c npc==3 
            for j in range(6): # 6 entries: ca, ha, va, cb, hb, vb
                pccomp[i][j]=stdcm[idx][i]*vecmag[idx][i][j] # for individual

            # v
            ax[0][i].bar(idx-.2,pccomp[i][2],width=0.35,color=ifname[idx][2],edgecolor='black',linewidth=0.45)
            ax[0][i].bar(idx+.2,pccomp[i][5],width=0.35,color=ifname[idx][2],edgecolor='black',linewidth=0.45,hatch='///')
            # h 
            ax[1][i].bar(idx-.2,pccomp[i][1],width=0.35,color=ifname[idx][2],edgecolor='black',linewidth=0.45)
            ax[1][i].bar(idx+.2,pccomp[i][4],width=0.35,color=ifname[idx][2],edgecolor='black',linewidth=0.45,hatch='///')

            for j in range(2):
                ax[j][i].set_ylim([0,2.5])

            # subtraction             
            rat[i][0]=pccomp[i][2]-pccomp[i][5] # va-vb 
            rat[i][1]=pccomp[i][1]-pccomp[i][4] # ha-hb 

            for j in range(2):
                axsub[j][i].bar(idx,rat[i][j],width=.6,color=ifname[idx][2],linewidth=0.45,
                                edgecolor='k') 
                
            axsub[0][i].set_ylim([-0.75,1.6])
            axsub[1][i].set_ylim([-1.2,1.6])
            
    # end loop

    xlabel=[]
    for idx in range(nfile):
        xlabel.append(ifname[idx][1])
    xlist=list(range(0,len(xlabel)))
            
    # both: plot settings
    plt.minorticks_on()                                            
    for j in range(2):
        for i in range(3): # b/c npc==3 
            ax[j][i].tick_params(axis="both",labelsize=12,length=7,width=2)
            axsub[j][i].tick_params(axis="both",labelsize=12,length=7,width=2)

            ax[0][i].set_title("PC"+str(i+1))
            axsub[0][i].set_title("PC"+str(i+1))

            ax[j][i].xaxis.set_minor_locator(MultipleLocator(1)) # to set the same 
            ax[j][i].xaxis.set_major_locator(MultipleLocator(1))
            ax[j][i].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax[j][i].yaxis.set_major_locator(MultipleLocator(1))
            ax[j][i].axes.xaxis.set_ticklabels([])
            if i>0: ax[j][i].axes.yaxis.set_ticklabels([])
            if j==1:
                ax[j][i].set_xticks(xlist)
                ax[j][i].set_xticklabels(xlabel,fontsize=12)
            
            axsub[j][i].xaxis.set_minor_locator(MultipleLocator(1)) # to set the same 
            axsub[j][i].xaxis.set_major_locator(MultipleLocator(1))
            axsub[j][i].yaxis.set_minor_locator(MultipleLocator(0.5))
            axsub[j][i].yaxis.set_major_locator(MultipleLocator(1))
            axsub[j][i].axes.xaxis.set_ticklabels([])                        
            axsub[j][i].axhline(0, c='black', ls='--')
            if i>0: axsub[j][i].axes.yaxis.set_ticklabels([])
            if j==1:
                axsub[j][i].set_xticks(xlist)
                axsub[j][i].set_xticklabels(xlabel,fontsize=12)

    ax[0][0].set_ylabel(r'Triad center of mass ($\AA$)',fontsize=12)
    ax[1][0].set_ylabel(r'Hinge ($\AA$)',fontsize=12)
    axsub[0][0].set_ylabel(r'V$\alpha$-V$\beta$ ($\AA$)',fontsize=12)
    axsub[1][0].set_ylabel(r'H$\alpha$-H$\beta$ ($\AA$)',fontsize=12)
            
    plt.show()
    quit()
