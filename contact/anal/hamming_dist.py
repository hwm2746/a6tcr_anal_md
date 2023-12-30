# plot hamming distance

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION 

pwd0='./out/hamming_'
fstart=25000 # frames to include [fstart,end] for histogram 

# ifname = location, plot label, style (color,marker}          
ifname=[ \
         ['WTlow','sandybrown',1.0],\
         ['WThigh','darkgoldenrod',0.9],\
]

######################################################################
def read_hamming(idx):
   """
    Read cordyn files
   """
   X1,Y1 = [],[]
   
   ifn=pwd0+ifname[idx][0]+'.dat'
   for line in open(ifn,'r'):
       if not line.startswith("#"): 
           values = [float(s) for s in line.split()]
           X1.append(values[0])
           Y1.append(values[2]) # hamming distance value

   return X1,Y1
       

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    nfile=len(ifname)

    # plot setup
    fig,ax= plt.subplots() # Hamming distance trajectory 
    hfig, hax = plt.subplots() # Histogram

    fig.suptitle('Hamming distance')
    hfig.suptitle('Hamming distance histogram')    
    
    X1,Y1 = [],[]
    X=[ [] for i in range(nfile)]
    Y=[ [] for i in range(nfile)]    
    for idx in range(nfile):
        X1,Y1= read_hamming(idx)
        X[idx]=np.flip(X1)                
        Y[idx]=np.flip(Y1)        
        
        hax.hist(Y1[fstart:len(Y1)],bins=10,
                 color=ifname[idx][1],edgecolor="black",density=True,
                 label=ifname[idx][0],alpha=ifname[idx][2])
        
    # Hamming distance trajectory 
    for idx in range(nfile):
        ax.plot(X[idx],Y[idx],label=ifname[idx][0],
                 color=ifname[idx][1], linewidth=0.08, alpha=1)
       
    # Hamming distance trajectory  
    ax.xaxis.set_minor_locator(MultipleLocator(200))
    ax.xaxis.set_major_locator(MultipleLocator(400))    
    ax.yaxis.set_minor_locator(MultipleLocator(2.5))
    ax.yaxis.set_major_locator(MultipleLocator(5))                
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.yaxis.set_tick_params(width=2,length=7)
    ax.xaxis.set_tick_params(width=2,length=7)
    ax.set_ylim([0,13.5])
    ax.set_xlabel('Time (ns)',fontsize=12)
    ax.set_ylabel('Hamming Dist.',fontsize=12)    

    # Histogram
    hax.xaxis.set_minor_locator(MultipleLocator(2.5))
    hax.xaxis.set_major_locator(MultipleLocator(5))
    hax.xaxis.set_tick_params(width=2,length=7)
    hax.yaxis.set_tick_params(width=2,length=7)
    hax.tick_params(axis='both', which='major', labelsize=12)
    hax.set_xlim([0,14])
    hax.set_xlabel('Hamming Dist.',fontsize=12)
    hax.set_ylabel('Histo.',fontsize=12)
    
    
    plt.show()

    quit()
