# force versus time, with occupancy colors

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION

pwd0='../data/'
pwd1='/fvocc_'
pwd2='.txt'

# range for occupancy colors
rval=[10,30,25,40] 

dname=['pmhc', 'intranoc'] 

ifname=[ \
         ['WT0','none'],\
         ['WTlow','none'],\
         ['WThigh','black'],\
]

#######################################################

######################################################################
def read_data(ifn):
   """
    Read cordyn files
   """
   t0, f0, occ0 = [], [], []
   
   for line in open(ifn,'r'):
       if not line.startswith("#"):
           values = [float(s) for s in line.split()]
           t0.append(values[0])  # time
           f0.append(values[1])  # force
           occ0.append(values[2])  # occupancy

   return t0, f0, occ0

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    nfile=len(ifname)

    fig,ax=plt.subplots(2,nfile,figsize=(10,5))
    fig.suptitle(r'Instantaneous force (pN) v. time' "\n" "Legends represent total occupancy",fontsize=12) 

    for idx in range(nfile):    
        for ip in range(len(dname)):
            ifn=pwd0+ifname[idx][0]+pwd1+dname[ip]+pwd2

            try:
                with open(ifn) as f:
                    print('File present')
            except FileNotFoundError:
                print('File is not present')
                continue # skip this reading, go to next dname 
                        
            time, force, occ = read_data(ifn)

            if ip==0:
               g0=ax[ip][idx].scatter(time,force,c=occ,cmap='plasma',
                                      vmin=rval[0],vmax=rval[1],
                                      s=80,edgecolor=ifname[idx][1])
            else: 
               g0=ax[ip][idx].scatter(time,force,c=occ,cmap='plasma',
                                      vmin=rval[2],vmax=rval[3],
                                      s=80,edgecolor=ifname[idx][1])
            
            ax[ip][idx].tick_params(axis="both",labelsize=12,length=7,width=2)
          
            ax[ip][idx].yaxis.set_major_locator(MultipleLocator(15))
            ax[ip][idx].xaxis.set_major_locator(MultipleLocator(500))    
         
            ax[ip][idx].set_xlim([175,1350])
            ax[ip][idx].set_ylim([-3,55])

            if (ip==0):
               ax[ip][idx].axes.xaxis.set_ticklabels([])
            else:
               ax[ip][idx].set_xlabel('Time (ns)',fontsize=12)
               
            if (idx==0):
               if (ip==0):
                  fig.colorbar(g0,ax=ax[ip][idx],ticks=[10,20,30])  
               if (ip==1):
                  fig.colorbar(g0,ax=ax[ip][idx],ticks=[25,30,35,40])

    ax[0][0].set_ylabel(r'TCR$\alpha\beta$-pMHC',fontsize=12)
    ax[1][0].set_ylabel(r'Intra-TCR$\alpha\beta$' "\n" r'excl. C$\alpha$-C$\beta$ contacts',fontsize=12)
                  
    plt.show()
    quit()                  
