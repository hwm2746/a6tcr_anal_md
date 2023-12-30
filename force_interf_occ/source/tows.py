#!/usr/bin/python
########################################################
#
#  tows.py: Tug-of-war sampling for a single constrained atom
#
#  Read an input file containing positions of constrained atoms during MD
#  simulation and carry out tug-of-war sampling.
#
#  usage: 
#  python tows.py -b begin -e end -s skip -v vmd_out -o oprefx
#  
#  vmd_out: vmd file (default: None)
#  begin: beginning frame number  (default: 0)
#  end  : 1 + [last frame number to include] (default: nfile)
#  skip : skip between frame number (default: 1)
#  oprefx :  output fle prefix
# 
#######################################################

import sys
import argparse
import struct

import numpy as np
from numpy import loadtxt, zeros, average,sqrt,linalg
import math

########################
# MODIFIABLE SECTION 

# The following from dyn_cont.inp
hk = 1.0 # spring const from CHARMM conns harm (kcal/(mol*A^2))

## 12
x0 =  93.8553237
x1 = -93.4270303

dt=0.02 # in ns

ifname='pos.dat' # from out_get_pos.dat, with grep, sed, gawk: see readme.txt

############
kT =  4.1419 # kBT in pN*nm

pos_ref = np.zeros((2,3))
pos_ref[0,0]= x0
pos_ref[1,0]= x1

#######################################################
# functions
#######################################################
def get_args():
    """
    Process the command line options.
    """

    parser = argparse.ArgumentParser(description=
                                     'Usage: python tows.py -b begin -e end -s skip -v vmd_out -o oprefx')
    parser.add_argument('-b', metavar='frm_ini', default=0)
    parser.add_argument('-e', metavar='frm_fin', default =-1)
    parser.add_argument('-s', metavar='skip',default=1)
    parser.add_argument('-v', metavar='vmd_out',default='')
    parser.add_argument('-o', metavar='oprefx')
    parser.add_argument('-w', metavar='twindow',type=int,default=-1)

    args = parser.parse_args() # args is a dictionary

    frm_ini= int(vars(args)['b']) 
    frm_fin = int(vars(args)['e']) 
    skip= int(vars(args)['s']) 
    vmd= vars(args)['v']
    oprefx= vars(args)['o']
    twindow= vars(args)['w']

    return frm_ini,frm_fin,skip,vmd,oprefx,twindow

def read_input(inFile):
    """
    Read atomic positions.
    """
    ff = open(inFile,'r')
    a=[]
    pos = [ [] for i in range(2) ]
    while 1:
        line = ff.readline()
        if not line:
            break
        if ( '#' not in line ) and ( len(line) > 1 ): # non-blank or comment
            linedata = line.split()            
            a.append(float(linedata[1]))
    for i in range(0,len(a),7):
        pos[0].append(a[i+1:i+4])
        pos[1].append(a[i+4:i+7])
    nframe=len(pos[0])

    return pos, nframe

def decimate(pos, frm_ini, frm_fin,skip):
    """ Extract sub arrays """

    pos_e=[ [] for j in range(2)]
    if frm_fin < 0 : frm_fin=len(pos[0])
    nf_e = (frm_fin - frm_ini)// skip # integer floor
    for j in range(2):
        pos_e[j] = [ pos[j][i] for i in range(frm_ini,frm_fin,skip) ]

    print('# nframe to use:', nf_e)
    return pos_e, nf_e
    

def tow(pos_e,nf_e):
    """
    Carry out tug-of-war sampling and get force. Ref structure is the first frame.
    """
    f_tow = np.zeros((2,4)) # force[:,4]: magnitude of force
    avg = np.zeros((2,3)) 
    var = np.zeros((2,3)) 
    cov = np.zeros((2,3)) #yz,zx,xy

    a=np.asarray(pos_e) # for indexing

    for i in range(2):
        for j in range(3):
            a[i,:,j] = a[i,:,j] - pos_ref[i,j]

    for i in range(2):
        for j in range(3):
            avg[i,j] = np.average(a[i,:,j])
        for j in range(3):
            r = np.std(a[i,:,j]);  var[i,j] = r * r
        for j in range(nf_e):
            cov[i,0] = cov[i,0] + a[i,j,1]*a[i,j,2]
            cov[i,1] = cov[i,1] + a[i,j,2]*a[i,j,0]
            cov[i,2] = cov[i,2] + a[i,j,0]*a[i,j,1]
        cov[i,0] = cov[i,0] / float(nf_e) - avg[i,1]*avg[i,2]
        cov[i,1] = cov[i,1] / float(nf_e) - avg[i,2]*avg[i,0]
        cov[i,2] = cov[i,2] / float(nf_e) - avg[i,0]*avg[i,1]
        
        for j in range(3):
            if j == 0:
                r = cov[i,1]*avg[i,2]/var[i,2] + cov[i,2]*avg[i,1]/var[i,1]
            elif j == 1:
                r = cov[i,2]*avg[i,0]/var[i,0] + cov[i,0]*avg[i,2]/var[i,2]
            else : # j == 2
                r = cov[i,0]*avg[i,1]/var[i,1] + cov[i,1]*avg[i,0]/var[i,0]
            f_tow[i,j] = 10*kT*(avg[i,j]-r)/var[i,j] # 10: force in pN

        fx = f_tow[i,0]; fy = f_tow[i,1]; fz = f_tow[i,2]
        f_tow[i,3] = np.sqrt( fx*fx + fy*fy + fz*fz )

    return f_tow

#######################################################
def getF(pos_e,nf_e):
    """
    Get force from spring const times displacement
    """
    f_spring = np.zeros((2,4)) # force[4]: magnitude of force
    r_avg = np.zeros((2,3))

    # 2*: for charmm convention. hk in pN/A
    hk0 = 2 * hk * 4184 * 0.1 / 6.022141 

    a=np.asarray(pos_e) # for indexing

    for i in range(2):
        for j in range(3):
            r_avg[i,j] = np.average(a[i,:,j])
            
    for i in range(2):
        for j in range(3):
            f_spring[i,j] = hk0 * (r_avg[i,j]- pos_ref[i,j])
        fx = f_spring[i,0]; fy = f_spring[i,1]; fz = f_spring[i,2]
        f_spring[i,3] = np.sqrt( fx*fx + fy*fy + fz*fz )

    return f_spring, r_avg
                
#######################################################
def windowF(pos_e,nf_e,frm_ini,twindow,oprefx):
    """
    Get force from spring const times displacement in intervals twindow
    """
    f_spring = np.zeros((2,4)) # force[4]: magnitude of force
    r_avg = np.zeros((2,3))

    # 2*: for charmm convention. hk in pN/A
    hk0 = 2 * hk * 4184 * 0.1 / 6.022141 

    a=np.asarray(pos_e) # for indexing;
    # accg: pos_e is truncated [frm_ini to frm_fin]; update for loop below 
    print("a",len(a[0]))
    nf0=len(a[0,:,0])
    print(nf0)

    ff = open(oprefx+'_interval.dat','w')
    ff.write('#  t        fx         fy         fz    f_spring\n')

    for frm in range(0,nf0,twindow): # accg
        #     for frm in range(frm_ini,nf0,twindow):
        for i in range(2):
            for j in range(3):
                r_avg[i,j] = np.average(a[i,frm:(frm+twindow),j])
                f_spring[i,j] = hk0 * (r_avg[i,j]- pos_ref[i,j])
            fx = f_spring[i,0]; fy = f_spring[i,1]; fz = f_spring[i,2]
            f_spring[i,3] = np.sqrt( fx*fx + fy*fy + fz*fz )

        fy0=0.5*(f_spring[0,1]-f_spring[1,1])
        fz0=0.5*(f_spring[0,2]-f_spring[1,2])
        fs0=0.5*(f_spring[0,3]+f_spring[1,3])
        t=dt* (float(frm+twindow/2) + frm_ini) #*dt # accg
        ff.write('{:7.2f} '.format(t))
        for i in range(3):
            f0=0.5*(f_spring[0,i]-f_spring[1,i])
            ff.write('{:9.2f} '.format(f0))
        ff.write('{:9.2f}\n'.format(fs0))

    return
                
#####################################################
# main routine
#####################################################
if __name__ == "__main__":
    frm_ini,frm_fin,skip,vmd,oprefx,twindow = get_args()

    pos,nframe = read_input(ifname)
    pos_e,nf_e = decimate(pos, frm_ini, frm_fin,skip)

    if twindow == -1:
        f_spring, r_avg = getF(pos_e,nf_e) # avg displacement * spring const
        f_tow = tow(pos_e,nf_e) # tug-of-war force at pos_e

        print('# Spring & TOWS force (pN):')
        for i in range(2):
            print('{0:+6.3e} {1:+6.3e} {2:+6.3e}  {3:+6.3e} {4:+6.3e} {5:+6.3e}'.format(f_spring[i,0],f_spring[i,1],f_spring[i,2],f_tow[i,0],f_tow[i,1],f_tow[i,2]))
            print('# Magnitude of Spring & TOWS force (pN):')
        for i in range(2):
            print('{0:6.3e} {1:6.3e}'.format(f_spring[i,3],f_tow[i,3]))

    else: # measure in windows
        windowF(pos_e,nf_e,frm_ini,twindow,oprefx)        
        
    exit()
