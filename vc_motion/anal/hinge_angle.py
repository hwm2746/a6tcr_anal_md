# Plot hinge angles 

import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import scipy.stats as stats
import math

#######################################################
### MODIFIABLE SECTION 

fstart=24999 # plot frames [fstart,end] 
#fstart=0 # for full trajectory plot, in SI
nbin0=10 # scale bins

pwd0='../data/'
pwd1='/meas_'
pwd2='.dat'

ifname=[ \
         ['Tab','purple','P',7,'dotted',2.5,'purple',0,'purple',0.1],\
         ['WT0','gold','v',6,'solid',5,'black',1.0,'gold',0.1],\
         ['WTlow','sandybrown','d',7,'solid',3.5,'black',1.5,'white',0.1],\
         ['WThigh','darkgoldenrod','o',5,'solid',3,'black',0.7,'darkgoldenrod'],\
]

#######################################################
def get_ref(): # first read all entries, get one w/ smallest angle distrib.

   dBin = [ 1000000 for i in range(2) ] # dBin per triad arm
   minsys = [ 'a' for i in range(2) ] # dBin per triad arm
   
   for idx in range(nfile):
      Y0,Y1=[],[]
      
      ifn=pwd0+ifname[idx][0]+pwd1+ifname[idx][0]+pwd2
      for line in open(ifn,'r'):
         if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            Y0.append(values[0])  # phi_a
            Y1.append(values[1])  # phi_b

      Y0=Y0[fstart:len(Y0)]
      Y1=Y1[fstart:len(Y1)]
         
      rdum=(max(Y0)-min(Y0))/nbin0
      #print("Y0 ",ifn, " rdum:",rdum)      
      if (rdum<dBin[0]):
         dBin[0]=rdum
         minsys[0]=ifn
         
      rdum=(max(Y1)-min(Y1))/nbin0
      #print("Y1 ",ifn, " rdum:",rdum)
      if (rdum<dBin[1]):
         dBin[1]=rdum
         minsys[1]=ifn
         

   for jp in range(2):
      print("minsys:",minsys[jp]," dBin:",dBin[jp])

   assert (dBin[0]!=1000000) 

######################################################################
def read_angle(idx):
   """
    Read cordyn files
   """
   Y0,Y1,Y2,Y3 = [], [], [], []

   ifn=pwd0+ifname[idx][0]+pwd1+ifname[idx][0]+pwd2
   
   for line in open(ifn,'r'):
       if not line.startswith("#"): 
           values = [float(s) for s in line.split()]
           Y0.append(values[0])  # phi_a
           Y1.append(values[1])  # phi_b
           Y2.append(values[5])  # psi_a
           Y3.append(values[6])  # psi_b

   return Y0, Y1, Y2, Y3

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
   nfile=len(ifname)
   #Y = [] * nfile
   Y0, Y1, Y2, Y3 = [] , [], [], []
   time=[ [] for i in range(nfile)]
   phi_a=[ [] for i in range(nfile)]
   phi_b=[ [] for i in range(nfile)]
   psi_a=[ [] for i in range(nfile)]
   psi_b=[ [] for i in range(nfile)]        

   ##############################
   # 1. execute this to get binning scheme per arm
   #get_ref()     
   # output:
   #minsys: /home/anacristina/Documents/md/tcr/1qse/nl/anal/boc/data_cab/meas_1qsenl.dat  dBin: 1.3635260000000002
   #minsys: /home/anacristina/Documents/md/tcr/1ao7/anal1/08/boc/data_cab/meas_1ao708.dat  dBin: 1.1137319999999988
   #
   # 2. now execute this w/ rest of code, after selecting a ref structure from ^ step 
   #dBin=get_bin_from_ref(2)
   #for i in range(2): print(dBin[i])
   dBin=[ -1 for i in range(2)]
   dBin[0]=1.3635260000000002
   dBin[1]=1.1137319999999988
         
   #############################   
   
   # plot variables
   fig,ax= plt.subplots(2) # 3 angle measures
   hfig,hax = plt.subplots(1,2,figsize=(13,5)) # histograms of above
   
   fig.suptitle(r'Hinge angles ($^\circ$)', fontsize=12)    
   ax[0].set_ylabel(r'$\angle$TCR$\alpha$',fontsize=12) 
   ax[1].set_ylabel(r'$\angle$TCR$\beta$',fontsize=12)

   hfig.suptitle(r'Histograms: Hinge angles ($^\circ$)', fontsize=12)
   hax[0].set_xlabel(r'$\angle$TCR$\alpha$',fontsize=12) 
   hax[1].set_xlabel(r'$\angle$TCR$\beta$',fontsize=12)
   
   for idx in range(nfile):
      Y0,Y1,Y2,Y3= read_angle(idx)
   
      Y0=Y0[fstart:len(Y0)]
      Y1=Y1[fstart:len(Y1)]
      Y2=Y2[fstart:len(Y2)]
      Y3=Y3[fstart:len(Y3)]
   
      ## print avgs
      #print(ifname[idx][1]," average va-ca:", stat.mean(Y0))
      #print(ifname[idx][1]," average vb-cb:", stat.mean(Y1))       

      ##########
      time0=list(range(0,len(Y0)))
      time = [ (x+fstart)/50 for x in time0]
      
      #for idx in range(nfile):
      ax[0].plot(time,Y0,label=ifname[idx][0],
                 color=ifname[idx][1], linewidth=0.1)
      
      ax[1].plot(time,Y1,label=ifname[idx][0],
                 color=ifname[idx][1], linewidth=0.1)

      ##########
      ## for va-ca 

      # calc nbin
      setrange=max(Y0)-min(Y0)
      nbin=math.ceil(setrange/dBin[0])
      #print("ifn:",idx," va-ca, bins=",nbin)

      ##
      density=stats.gaussian_kde(Y0)
      n,x, _ = hax[0].hist(Y0,bins=50,
                            histtype=u'step',density=True,
                           color=ifname[idx][1],alpha=0)
      hax[0].plot(x,density(x),color=ifname[idx][1],
                  linestyle=ifname[idx][4],linewidth=ifname[idx][5])
   
      data = np.histogram(Y0,bins=nbin,density=True)
      ## data=[hist value, bin value, patches]
      for ibin in range(nbin): 
         xval=(data[1][ibin]+data[1][ibin+1])/2.0
         hax[0].plot(xval,data[0][ibin],color=ifname[idx][1],
                     marker=ifname[idx][2],markerfacecolor=ifname[idx][8],
                     markersize=ifname[idx][3]*1.2, markeredgewidth=ifname[idx][7],
                     markeredgecolor=ifname[idx][6])         

      #####
      ## for vb-cb 

      # calc nbin
      setrange=max(Y1)-min(Y1)
      nbin=math.ceil(setrange/dBin[1])
      #print("ifn:",idx," vb-cb, bins=",nbin)

      ##
      density=stats.gaussian_kde(Y1)
      n,x, _ = hax[1].hist(Y1,bins=50,
                            histtype=u'step',density=True,
                            color=ifname[idx][1],alpha=0)
      hax[1].plot(x,density(x),color=ifname[idx][1],
                  linestyle=ifname[idx][4],linewidth=ifname[idx][5])
   
      data = np.histogram(Y1,bins=nbin,density=True)
      ## data=[hist value, bin value, patches]
      for ibin in range(nbin): 
         xval=(data[1][ibin]+data[1][ibin+1])/2.0
         hax[1].plot(xval,data[0][ibin],color=ifname[idx][1],
                     marker=ifname[idx][2],markerfacecolor=ifname[idx][8],
                     markersize=ifname[idx][3]*1.2, markeredgewidth=ifname[idx][7],
                     markeredgecolor=ifname[idx][6])        
         
      plt.minorticks_on()
      for jp in range(2):       
   
         ax[jp].tick_params(axis="both",labelsize=12,length=7,width=2)

         ax[jp].xaxis.set_minor_locator(MultipleLocator(100))
         ax[jp].xaxis.set_major_locator(MultipleLocator(200))  
         ax[jp].yaxis.set_major_locator(MultipleLocator(10))
         ax[jp].yaxis.set_minor_locator(MultipleLocator(5))
   
         hax[jp].xaxis.set_major_locator(MultipleLocator(10))
         hax[jp].xaxis.set_minor_locator(MultipleLocator(5))
         hax[jp].yaxis.set_major_locator(MultipleLocator(0.1))
         hax[jp].yaxis.set_minor_locator(MultipleLocator(0.05))
          
         hax[jp].tick_params(axis="both",labelsize=12,length=7,width=2)
   
   ax[0].set_ylim([58,88])
   ax[1].set_ylim([64,83])
   
   hax[0].set_xlim([58,88])
   hax[1].set_xlim([64,83])
   hax[0].set_ylim([-0.005,0.3])
   hax[1].set_ylim([-0.005,0.3])
   
   ax[1].set_xlabel('Time (ns)',fontsize=12)         
   hax[0].set_ylabel('Norm. Distrib.',fontsize=12)
   
   plt.show()
   quit()
