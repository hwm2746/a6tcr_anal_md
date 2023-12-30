# Read output pca_dotp.cpp
# Plot color map of dot prod values between systems 

import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import networkx as nx

#######################################################
### MODIFIABLE SECTION 

ip0=[0,1,2]

names=[r'V$\alpha\beta$',r'T$\alpha\beta$','V','0','Low','High'] # order should match setSys mols list in pca_dotp.cpp 

fname='out_dotp.txt'

#######################################################

######################################################################
def read_data(ifn,pccf,ip):
    """
    Read pcprod data
    """
    
    for line in open(ifn,'r'):
        if not line.startswith("#"):
            values= [float(s) for s in line.split()]
            i=int(values[0])
            j=int(values[1])
            pcmod=int(values[2])
            if (pcmod==ip):
                pccf[j][i]=abs(values[3]) # for viewing purposes
                
    return pccf

#######################################################
# Main routine
#######################################################

if __name__ == "__main__":
    
    pccf = [[100]* len(names) for i in names]
    fig,ax=plt.subplots(1,3,figsize=(15,5))
    fig.suptitle('Triad PC dot product', fontsize=12)

    for ip in range(len(ip0)): 
        pccf=read_data(fname,pccf,ip)
    
        ########### plot ##########
        g0 = img = ax[ip].imshow(pccf,'plasma',vmax=1.0,vmin=0.0) 
    
        for (j,i),label in np.ndenumerate(pccf):
            if (label>1.1): continue        
            if (float(label)>0.5): tcol="black"
            else : tcol="white"
            ax[ip].text(i,j,"%.2f"%round(label,2),ha='center',va='center',
                           fontsize=12,color=tcol)
            
    # end     for ip in range(len(ip0)): 
            
    xlist = list(range(0,len(names)))            

    for ip in range(len(ip0)):
        ax[ip].title.set_text('PC'+str(ip+1))
        ax[ip].set_xticks(xlist)
        ax[ip].set_xticklabels(names,fontsize=12)
        ax[ip].set_yticks(xlist)
        ax[ip].set_yticklabels(names,fontsize=12)

    
    plt.show()
    quit()
