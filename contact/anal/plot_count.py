# Plot sum of TCRab-pMHC contact counts

import math 
import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#######################################################
### MODIFIABLE SECTION 

pwd0='./out/ncont_'
pwd1='.dat'

dname=['ma','mb','pa','pb']

# color of bars. follow order in ./out/ncont_*.dat. 
pmhc_scol=['gold','sandybrown','darkgoldenrod','red']               


#######################################################

#######################################################
def read_data(ifn):
    
    sys,Yavg,Ystd=[],[],[]    
    
    for line in open(ifn,'r'):
        if not line.startswith("#"):
            list_of_strings=line.split()
            sys.append(list_of_strings[0])
            values = [float(s) for s in list_of_strings[1:]]
            Yavg.append(values[0])
            Ystd.append(values[1])
            
    return sys,Yavg, Ystd 

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":

    sys=[ [] for i in range(len(dname))]    
    Y=[ [] for i in range(len(dname))]
    std=[ [] for i in range(len(dname))]

    for ip in range(len(dname)):
        ifn=pwd0+dname[ip]+pwd1
        try:
            with open(ifn) as f:
                print('File present')
        except FileNotFoundError:
            print('File is not present')
            continue # skip this reading, go to next dname 
        
        sys0,Y0,std0 = read_data(ifn)

        for j in range(len(sys0)): # save individually
            sys[ip].append(sys0[j])
            Y[ip].append(Y0[j])
            std[ip].append(std0[j])

    ## sum interface contacts
    pmhc = [0] * len(sys[0]) # systems w/ pmhc
    var_pmhc = [0] * len(sys[0]) # systems w/ pmhc
    std_pmhc = [0] *len(sys[0]) # average std

    for ip in range(len(dname)):
        if (dname[ip]=='pa' or dname[ip]=='pb'
            or
            dname[ip]=='ma' or dname[ip]=='mb'):

            for j in range(len(sys[ip])):
                pmhc[j] = Y[ip][j] + pmhc[j]
                var_pmhc[j] = (std[ip][j]*std[ip][j]) + var_pmhc[j]  # pmhc variances

    for j in range(len(sys[0])): # average variances
        std_pmhc[j] = math.sqrt(var_pmhc[j] / len(sys[0])) # 4 variances summed 

    #print('pmhcs:',std_pmhc)
    # for std: average the variances; then take square root
    # to get the average standard deviation. 
    
    ##############################
    # PLOT
    
    fig = plt.figure()
    ax=plt.axes()
    plt.title('Sum of contacts with overall occ > 50%, inst. occ > 80%') 

    bars = ax.bar(sys[0],pmhc,color=pmhc_scol,edgecolor='k',linewidth=1.5,width=.6)

    ax.errorbar(sys[0],pmhc,yerr=std_pmhc,linestyle='',
                solid_capstyle='projecting',capsize=5,color='black',
                elinewidth=2,capthick=2) 
    plt.setp(ax.get_xticklabels()) # show x labels

    plt.minorticks_on()
    ax.tick_params(axis="y",labelsize=12,length=10,width=3)
    ax.tick_params(axis="x",labelsize=12,length=10,width=3)
    ax.yaxis.set_major_locator(MultipleLocator(5.0))
    ax.yaxis.set_minor_locator(MultipleLocator(2.5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))        

    ax.set_ylim([7,28])
    ax.set_ylabel('Num. contacts with pMHC',fontsize=12)

    plt.show()
    quit()
                

             
