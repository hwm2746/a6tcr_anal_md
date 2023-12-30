""" occ_traj.py: Plot occupancy trajectory 
"""

import sys
import argparse
import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
import matplotlib.ticker as ticker # for scaling tick values

############
ofmt='pdf' # output format (pdf, png, etc)
############


#######################################################
def get_args():
    """
    Command line options.

    suffix: suffix for occupancy file, in case occupancy is for an interval.
         ex) In data_*/[hb,np]_*0_occ.dat, suffix is 0
         To consider bond occupancy for the whole traj, omit suffix

    namedef: name definitions ("namedef_*.py" without ".py")
    ofname: output image filename without extension (pdf)
    """
    parser = argparse.ArgumentParser(description=
             'Usage: occ_traj. -i ifname -s suffix -o ofname -r [resi list]')
    #parser.add_argument('-i', metavar='ifname')
    parser.add_argument('-s', metavar='suffix')
    parser.add_argument('-n', metavar='namedef')
    parser.add_argument('-o', metavar='ofname')

    args = parser.parse_args() # args is a dictionary
    #ifname = vars(args)['i']
    namedef = vars(args)['n']
    ofname = vars(args)['o']
    suffix = vars(args)['s']
    if suffix == None: suffix=''
    return namedef, ofname, suffix

#######################################################
def three2one(resn):
    """  
    convert 3-letter to 1-letter code of amino acid
    """
    
    sdum='XXX'
    if (resn=='ALA'): sdum='A'
    elif (resn=='ARG'): sdum='R'
    elif (resn=='ASN'): sdum='N'
    elif (resn=='ASP'): sdum='D'
    elif (resn=='CYS'): sdum='C'
    elif (resn=='GLN'): sdum='Q'
    elif (resn=='GLU'): sdum='E'
    elif (resn=='GLY'): sdum='G'
    elif (resn=='HIS'): sdum='H'
    elif (resn=='HSD'): sdum='H'
    elif (resn=='HSE'): sdum='H'
    elif (resn=='ILE'): sdum='I'
    elif (resn=='LEU'): sdum='L'
    elif (resn=='LYS'): sdum='K'
    elif (resn=='MET'): sdum='M'
    elif (resn=='PHE'): sdum='F'
    elif (resn=='PRO'): sdum='P'
    elif (resn=='SER'): sdum='S'
    elif (resn=='THR'): sdum='T'
    elif (resn=='TRP'): sdum='W'
    elif (resn=='TYR'): sdum='Y'
    elif (resn=='VAL'): sdum='V' 
    elif (resn=='HYP'): sdum='O'
    else:
        print('Unrecognized argument: ',resn)

    return sdum

#######################################################
def read_bond(suffix,ifname,nfile):
    """
    Read occupancy-ranked bond data.
    """
    occ_hb=[ [] for i in range(nfile)]
    occ_np=[ [] for i in range(nfile)]

    for i in range(nfile):
        ifn0=ifpwd0+'/data_'+ifname[i]+'/'
        ifn1=ifname[i]+suffix+'_occ.dat'
        ifn=ifn0+'hb_'+ifn1 # hbond
        with open(ifn, 'r') as ff:
            for line in ff:
                if line.startswith('#'): continue # skip comment
                occ_hb[i].append(line.split())
        ff.close()
        ifn=ifn0+'np_'+ifn1 # nonpolar
        with open(ifn, 'r') as ff:  
            for line in ff:
                if line.startswith('#'): continue # skip comment
                occ_np[i].append(line.split())
        ff.close()

    return occ_hb,occ_np

#######################################################
def sele_bond1(occ_hb,occ_np,ocut,ocut1,ifname,nfile):
    """ 
    Select bond to process, with occupancy during any time of the 
    trajectory higher than ocut AND overall occupancy > ocut1
    """
    sele_hb=[ [] for i in range(nfile)]
    sele_np=[ [] for i in range(nfile)]

    for i in range(nfile):
        ifn0=ifpwd0+'/data_'+ifname[i]+'/'        
        ifn1=ifname[i]+'_occ_traj.dat'
        ifn_hb=ifn0+'hb_'+ifn1; ifn_np=ifn0+'np_'+ifn1
        if (len(occ_hb[i])>0): # process hbond
            a0=np.loadtxt(ifn_hb)
            if len(occ_hb[i]) != len(a0[0])-1: # a0[0] has frame number col.
                print('ERROR: H-bond number mismatch.')
                quit()
            # print 'Num. hbond for ',ifname[i]
            for j in range(1,len(a0[0])): # go over hb
                if (float(occ_hb[i][j-1][-1]) > ocut1):            
                    q=int(occ_hb[i][j-1][0])
                    for k in range(len(a0)): # go over frames
                        #if a0[k][j] >ocut :
                        if a0[k][q+1] >ocut :
                            # j-1: account for frame number col
                            sele_hb[i].append(int(occ_hb[i][j-1][0]))
                            #print(ifn_hb,j,q) # for debug
                            break
        else:
            print('No hbonds for ',ifname[i])

        if (len(occ_np[i])>0): # process nonpolar
            a1=np.loadtxt(ifn_np)
            if len(occ_np[i]) != len(a1[0])-1:  # a1[0] has frame number col.
                print('ERROR: Nonpolar bond number mismatch.')
                quit()
            for j in range(1,len(a1[0])): # go over np
                if (float(occ_np[i][j-1][-1]) > ocut1):
                    q=int(occ_np[i][j-1][0])
                    for k in range(len(a1)): # go over frames
                        #if a1[k][j] >ocut:
                        if a1[k][q+1] >ocut:
                            sele_np[i].append(int(occ_np[i][j-1][0]))
                            
                            break
        else:
            print('No nonpolar contacts for ',ifname[i])
        #print ifn0, a0.shape, len(a0)
    # quit()
    
    return sele_hb,sele_np

#######################################################
def get_label(sele_hb,sele_np,occ_hb,occ_np):
    """ 
    Prepare bond labels
    """
    label_hb=[ [] for i in range(nfile)]
    label_np=[ [] for i in range(nfile)]

    for i in range(nfile):
        if len(sele_hb[i]) == 0 : continue
        for j in sele_hb[i]:
            sdum=''
            for k in range( len(occ_hb[i]) ):
                if int(occ_hb[i][k][0]) == j:
                    sdum=three2one( occ_hb[i][k][2] )+str(int(occ_hb[i][k][3]))\
                    +'-' +three2one( occ_hb[i][k][5] )+str(int(occ_hb[i][k][6]))
                    break
            if len(sdum) == 0:
                print(ifname[i],': no label found for bond ',j)
            label_hb[i].append(sdum)

    for i in range(nfile):
        if len(sele_np[i]) == 0 : continue
        for j in sele_np[i]:
            sdum=''
            for k in range( len(occ_np[i]) ):
                if int(occ_np[i][k][0]) == j:
                    sdum=three2one( occ_np[i][k][2] )+str(int(occ_np[i][k][3]))\
                    +'-' +three2one( occ_np[i][k][5] )+str(int(occ_np[i][k][6]))
                    break
            if len(sdum) == 0:
                print(ifname[i],': no label found for bond ',j)
            label_np[i].append(sdum)

    return label_hb, label_np

#######################################################
def plot_occ(sele_hb,sele_np,occ_hb,occ_np,ifname,
             ofname,gname,molname,nfile):
    """ 
    Plot occupancy trajectory.
    """

    # https://matplotlib.org/tutorials/intermediate/constrainedlayout_guide.html#sphx-glr-tutorials-intermediate-constrainedlayout-guide-py
    pc_kwargs = {'rasterized': True, 'cmap': 'Spectral'} # , 'norm': norm}
    # viridis, Spectral

    label_hb, label_np = get_label(sele_hb,sele_np,occ_hb,occ_np)
    #molname=ofname.upper()
    
    ############
    # h-bond
    ############

    nhb=[]; nplot=0
    for i in range(nfile): # count number of bond files to process
        j=len(sele_hb[i])
        if j == 0 : continue # skip entry with no bond count
        nhb.append(j)
        nplot=nplot+1

    if not (len(nhb)==0):
        #nhb[0]=len(nhb)*0.2
        gs = gridspec.GridSpec(1, 1) #, height_ratios=nhb)
        w0=10; #h0=0.3*sum(nhb)+1/sum(nhb)
        h0=0.5*sum(nhb)
        #w0=10; h0=0.3*10+0.7
        fig = plt.figure(figsize=(w0,h0)) # original
        g_id=0 # graph counter
        for i in range(nfile):
            if len(sele_hb[i]) == 0 : continue
            ifn0=ifpwd0+'/data_'+ifname[i]+'/'        
            ifn1=ifname[i]+'_occ_traj.dat'
            ifn=ifn0+'hb_'+ifn1
            a0=np.loadtxt(ifn)
            # a1: increment column number by 1 to account for time col in occ
            a1=a0[:,[ x + 1 for x in sele_hb[i]] ]
            nframe=len(a1)
            n=len(a1[0,:])
            a2=np.transpose(a1)
            ax = fig.add_subplot(gs[g_id])
            if g_id == 0 : # top graph
                ax.set_title(molname+' (H-bond)',fontweight='bold')
            ax.set_xticks(np.arange(0,nframe,10000) )
            # bond label
            plt.yticks(np.arange(0.5,len(label_hb[i])+.5,1),label_hb[i])
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.set_ylabel(gname[i],size=10,rotation=90,fontweight='bold')
            # get y offset: since height is 0.3*(len(sele_hb[i]) inches,
            # 10pt font =13pixels=0.18in, subtract 0.28/height from 0.5
            # (0.28: chosen by try/error)
            ys0=0.5-0.28/(0.3*float(len(sele_hb[i])))
            ax.yaxis.set_label_coords(-0.18,ys0)
            #ax.vlines(-0.15,0,1)#,'k-',0.02)
            #xs0=-0.15-float(g_id%2)*0.03 # zig-zag bond labels
            #ax.yaxis.set_label_coords(xs0,0.5)
            
            # xtick label only for the bottom graph
            if g_id < (nplot-1): ax.set_xticklabels('')
            else:
                # change from steps to ns (0.02ns= nsavc*step)
                # https://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
                xnew = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.02))
                ax.xaxis.set_major_formatter(xnew)
                
                ax.set_xlabel('Time (ns)',size=10,fontweight='bold')
                
            im = ax.pcolormesh(a2, **pc_kwargs)
            g_id = g_id + 1; # graph counter

        fig.subplots_adjust(\
                            left  = 0.2  ,  # left margin
                            right =  0.8  ,  # start of right margin
                            bottom =0.35 ,  # bottom margin
                            top =  0.8    ,  # start of top margin
                            wspace = 0.05 ,  # hgap
                            hspace = w0*0.05/h0   # vgap: make prop to aspect ratio
                            )
        
        # a single colorbar
        # https://stackoverflow.com/questions/41860386/one-colorbar-for-multiple-pandas-subplots
        
        ## Uncomment to show colorbar
        #im = plt.gca().get_children()[0]
        #cbar = fig.add_axes([0.82,0.4,0.015,0.3])  #x,y,w,h
        #cbar.tick_params(labelsize=10)
        #fig.colorbar(im, cax=cbar) #, ticks=np.arange(0,1,.25))
        #cbar.set_ylabel('Occupancy',size=10) # this should be below fig.colorbar
        
        fig.savefig(ofname+'_hb_occ.'+ofmt, format=ofmt)

    ############
    # nonpolar
    ############
    
    nnp=[]; nplot=0    
    for i in range(nfile): # count number of bond files to process
        j=len(sele_np[i])
        if j == 0 : continue # skip entry with no bond count
        nnp.append(j)
        nplot=nplot+1

    #nnp[0]=0.32
    if not (len(nnp)==0):
        gs = gridspec.GridSpec(1, 1)#, height_ratios=nnp)
        #w0=10; h0=0.3*sum(nnp)+1
        w0=10; #h0=0.3*sum(nnp)+1/sum(nnp)
        h0=0.5*sum(nnp)        
        fig = plt.figure(figsize=(w0,h0 )) # original
        #fig = plt.figure(figsize=(10,3)) # hard-set 

        g_id=0 # graph counter
        for i in range(nfile):
            if len(sele_np[i]) == 0 : continue


            ifn0=ifpwd0+'/data_'+ifname[i]+'/'
            ifn1=ifname[i]+'_occ_traj.dat'
            ifn=ifn0+'np_'+ifn1
            a0=np.loadtxt(ifn)
            # a1: increment column number by 1 to account for time col in occ
            a1=a0[:,[ x + 1 for x in sele_np[i]] ]
            nframe=len(a1)
            n=len(a1[0,:])
            a2=np.transpose(a1)
            ax = fig.add_subplot(gs[g_id])
            if g_id == 0 : # top graph
                ax.set_title(molname+' (Nonpolar)',fontweight='bold')
            ax.set_xticks(np.arange(0,nframe,10000) )
            # bond label
            plt.yticks(np.arange(0.5,len(label_np[i])+.5,1),label_np[i])
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.set_ylabel(gname[i],size=10,rotation=90,fontweight='bold')
            ys0=0.5-0.28/(0.3*float(len(sele_np[i]))) # see H-bond above
            ax.yaxis.set_label_coords(-0.18,ys0)
        
        
            # xtick label only for the bottom graph
            if g_id < (nplot-1): ax.set_xticklabels('')
            else:
                # change from steps to ns (0.02ns= nsavc*step)
                # https://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
                xnew = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.02))
                ax.xaxis.set_major_formatter(xnew)
                
                ax.set_xlabel('Time (ns)',size=10,fontweight='bold')
        
            im = ax.pcolormesh(a2, **pc_kwargs)
            g_id = g_id + 1; # graph counter

        fig.subplots_adjust(\
                            left  = 0.2  ,  # left margin
                            right = 0.8  ,  # start of right margin
                            bottom = 0.35 ,  # bottom margin
                            top = 0.85    ,  # start of top margin
                            wspace = 0.05 ,  # hgap
                            hspace = w0*0.05/h0   # vgap: make prop to aspect ratio
                            # hspace = 0.1   # vgap
                            )

        # a single colorbar
        # https://stackoverflow.com/questions/41860386/one-colorbar-for-multiple-pandas-subplots
        
        ## Uncomment to show colorbar
        #im = plt.gca().get_children()[0]
        #cbar = fig.add_axes([0.82,0.4,0.015,0.3])  #x,y,w,h
        #cbar.tick_params(labelsize=10)
        #fig.colorbar(im, cax=cbar) #, ticks=np.arange(0,1,.25))
        #cbar.set_ylabel('Occupancy',size=10) # this should be below fig.colorbar
    
        fig.savefig(ofname+'_np_occ.'+ofmt, format=ofmt)
        
        ### END NONPOLAR ### 

    quit()

    

#######################################################
# Main routine
#######################################################
if __name__ == "__main__":
    namedef,ofname,suffix = get_args()

    nd= __import__(namedef) # from a import *
    nfile=len(nd.ifname)
    ifn0=nd.ifname
    gname0=nd.gname
    title0=nd.molname
    ifpwd0=nd.ifpwd
    
    occ_hb,occ_np=read_bond(suffix,ifn0,nfile)
    sele_hb,sele_np=sele_bond1(occ_hb,occ_np,nd.ocut,nd.ocut1,ifn0,nfile)
    plot_occ(sele_hb,sele_np,occ_hb,occ_np,ifn0,ofname,gname0,title0,nfile)
    
    quit()
