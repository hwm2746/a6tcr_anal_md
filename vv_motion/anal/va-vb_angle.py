# Plot va-vb triad angles

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
pwd1='/phi'
pwd2='.txt'
triarm=['1','2','3'] # plot angles for all 3 triad arms (e1-e1, etc)
ifname=[ \
         ['Vab','blue','P',7,'solid',1.5,'blue',0,'blue'],\
         ['Tab','purple','P',7,'dotted',1.5,'purple',0,'purple'],\
         ['Vab-pMHC','red','o',6,'dashed',2,'k',0.3,'red'],\
         ['WT0','gold','v',5,'solid',4.5,'k',1.0,'gold'],\
         ['WTlow','sandybrown','d',8.5,'solid',4,'k',1.6,'white'],\
         ['WThigh','darkgoldenrod','o',7,'solid',5,'k',0.5,'darkgoldenrod'],\
]
#######################################################

#######################################################
def get_ref(): # first read all entries, get one w/ smallest angle distrib.

   dBin = [ 1000000 for i in range(3) ] # dBin per triad arm
   minsys = [ 'a' for i in range(3) ] # dBin per triad arm
   
   for idx in range(nfile):
      for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
         X0=[]
         ifn=pwd0+ifname[idx][0]+'/'+pwd1+triarm[jp]+pwd2
         for line in open(ifn,'r'):
            if not line.startswith("#"): 
               values = [float(s) for s in line.split()]
               X0.append(values[0])

         X0=X0[fstart:len(X0)]
         rdum=(max(X0)-min(X0))/nbin0
         if (rdum<dBin[jp]):
            dBin[jp]=rdum
            minsys[jp]=ifn

   assert (dBin[0]!=1000000) 

######################################################################
def get_bin_from_ref(idx): # first read all entries, get one w/ smallest angle distrib.

   dBin = [ 1000000 for i in range(3) ] # dBin per triad arm
   
   for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
      X0=[]
      ifn=pwd0+ifname[idx][0]+'/'+pwd1+triarm[jp]+pwd2
      for line in open(ifn,'r'):
         if not line.startswith("#"): 
            values = [float(s) for s in line.split()]
            X0.append(values[0])
            
      X0=X0[fstart:len(X0)]
      dBin[jp]=(max(X0)-min(X0))/nbin0

   assert (dBin[0]!=1000000) 
   return dBin

######################################################################
def read_angle(idx,jp):
   """
    Read triad arm angle files
   """
   Y0=[]

   # triad arm angles
   ifn=pwd0+ifname[idx][0]+'/'+pwd1+triarm[jp]+pwd2
   #print(ifn)
   
   for line in open(ifn,'r'):
      if not line.startswith("#"): 
         values = [float(s) for s in line.split()]
         Y0.append(values[0])
         
   return Y0

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":

   nfile=len(ifname)
   time0=[ [] for i in range(nfile)]
   time=[ [] for i in range(nfile)]      
   Y=[ [] for i in range(nfile)]     # e3 triad arm angle

   #############################
   # 1. execute this to get binning scheme per arm
   #get_ref()     
   # output:
   #minsys: 1ao7/pmv/anal/boc/phi1.txt  dBin: 1.8091999999999984
   #minsys: 1ao7/pmv/anal/boc/phi2.txt  dBin: 1.6252899999999997
   #minsys: 1ao7/pmv/anal/boc/phi3.txt  dBin: 1.7116100000000003
   
   ### 2. now execute this w/ rest of code, after selecting a ref structure from ^ step 
   ##dBin=get_bin_from_ref(2)
   dBin=[ -1 for i in range(3)]
   dBin[0]=1.8091999999999984
   dBin[1]=1.6252899999999997
   dBin[2]=1.7116100000000003
   
   ##############################   
   
   # subplots(nfiles, npc modes)
   fig,ax = plt.subplots(3,figsize=(10,5)) # 3 angle measures
   fig.suptitle(r'V$\alpha$-V$\beta$ triad angles', fontsize=12)

   hfig,hax = plt.subplots(1,3,figsize=(15,5)) # histograms
   hfig.suptitle(r'Histograms: V$\alpha$-V$\beta$ triad angles', fontsize=12)   
   
   for idx in range(nfile): 
      for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
         
         ax[jp].set_ylabel(r'$\angle$e'+str(jp+1)+' ($^\circ$)',fontsize=12)
         hax[jp].set_xlabel(r'$\angle$e'+str(jp+1)+' ($^\circ$)',fontsize=12)
   
         Y1=[]
         Y1=read_angle(idx,jp)
         Y[idx]=Y1[fstart:len(Y1)]
         time0[idx]=list(range(0,len(Y[idx])))
         time[idx] = [ (x+fstart)/50 for x in time0[idx]]
         #myList[:] = [x / myInt for x in myList]
   
         ## print avg,std
         #print(ifname[idx][1]," average angle:", stat.mean(Y[idx]))
         #print(ifname[idx][1]," std:", np.std(Y[idx]))       
   
         
         ax[jp].plot(time[idx],Y[idx],#label=ifname[idx][1],
                     color=ifname[idx][1],linewidth=0.1) # ifname[idx][10])
   
         #####
         # calc nbin
         setrange=max(Y[idx])-min(Y[idx])
         nbin=math.ceil(setrange/dBin[jp])
         #print("ifn:",idx," arm:",jp," bins=",nbin)
         #####
         
         density=stats.gaussian_kde(Y[idx])
         n,x, _ = hax[jp].hist(Y[idx],bins=50,
                               histtype=u'step',density=True,
                               color=ifname[idx][1],alpha=0)
         hax[jp].plot(x,density(x),color=ifname[idx][1],
                      linestyle=ifname[idx][4],linewidth=ifname[idx][5])
   
         data = np.histogram(Y[idx],bins=nbin,density=True)
         ## data=[hist value, bin value, patches]
         for ibin in range(nbin): 
            xval=(data[1][ibin]+data[1][ibin+1])/2.0
            hax[jp].plot(xval,data[0][ibin],color=ifname[idx][1],
                         marker=ifname[idx][2],markerfacecolor=ifname[idx][8],
                         markersize=ifname[idx][3], markeredgewidth=ifname[idx][7],
                         markeredgecolor=ifname[idx][6])
         plt.minorticks_on()
         
         ax[jp].xaxis.set_minor_locator(MultipleLocator(100))
         ax[jp].xaxis.set_major_locator(MultipleLocator(200))         
         ax[jp].yaxis.set_minor_locator(MultipleLocator(5))
         ax[jp].yaxis.set_major_locator(MultipleLocator(10))
         ax[jp].tick_params(axis="both",labelsize=12,length=7,width=2)
         if (jp<2): 
            ax[jp].axes.xaxis.set_ticklabels([])
            
   ax[2].set_xlabel('Time (ns)',fontsize=12)         
   hax[0].set_ylabel('Norm. Distrib.',fontsize=12)
 
   #### end loop
   
   for jp in range(len(triarm)):  # triad arm 1, 2, 3.. etc
      hax[jp].xaxis.set_minor_locator(MultipleLocator(7.5))
      hax[jp].xaxis.set_major_locator(MultipleLocator(15))
      
      #hax[jp].set_xticks([])
      if (jp>0):
         hax[jp].set_yticklabels([])
         hax[jp].yaxis.set_minor_locator(MultipleLocator(.05))
         hax[jp].yaxis.set_major_locator(MultipleLocator(0.1))
      hax[jp].tick_params(axis="both",labelsize=12,length=7,width=2) 
   
   ## WT
   hax[0].set_xlim([143,181]) 
   hax[1].set_xlim([76,111])
   hax[2].set_xlim([70,103])
   hax[0].set_ylim([-0.005,0.23])
   hax[1].set_ylim([-0.005,0.23])
   hax[2].set_ylim([-0.005,0.23])
   
   # for angle trajectory, use same as hax above
   ax[0].set_ylim([143,181])
   ax[1].set_ylim([76,111])
   ax[2].set_ylim([70,105])
   
   plt.show()
   quit()

