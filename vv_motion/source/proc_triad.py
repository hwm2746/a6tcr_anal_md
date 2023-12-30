"""
proc_triad.py: process triad & build BOC
"""
import sys
import numpy as np

from name_def import get_names
segi, vn,cn,hn,molname, fdir =get_names()

def read_raw(ifname):
    """ read raw data """
    a=np.loadtxt(ifname,dtype=float)
    rd=[]; rc=[]
    for i in range(0,len(a),7):
        k=int(a[i])
        rd.append([a[i+1],a[i+2],a[i+3]])
        rc.append([a[i+4],a[i+5],a[i+6]])
    return rd, rc

def read_cm(ifname):
    """ read cm data """
    a=np.loadtxt(ifname,dtype=float)
    rc=[]
    for i in range(0,len(a),4):
        k=int(a[i])
        rc.append([a[i+1],a[i+2],a[i+3]])
    return rc

def triad_gen(rdomain):
    ifn0=fdir+'raw_'+rdomain+'0.dat'
    ifn1=fdir+'raw_'+rdomain+'1.dat'

    rda0, rca0=read_raw(ifn0)
    rda1, rca1=read_raw(ifn1)
    e0=0.5*np.add(rca0,rca1) # centroid

    e3=0.5*np.add(rda0,rda1)
    e3=normalize(e3)

    # build e2 & e3
    e10=np.subtract(rca0,rca1) # rca1-rca0: point outward 
    e10=normalize(e10)
    e2=np.cross(e3,e10); e2=normalize(e2)
    e1=np.cross(e2,e3); e1=normalize(e1)
    
    return [e0,e1,e2,e3]
    
def normalize(a):
    """ normalize array of vectors """
    for i in range(len(a)):
        n=np.linalg.norm(a[i])
        a[i]=a[i]/n
    return a

def write_psf(segname,dname):
    """ Generate ofname.psf containing triad info """
    ff = open(fdir+dname+".psf",'w')
    ntriad=1 # len(names) # number of triads to build
    natom = ntriad*4
    nbond = ntriad* 3
    ff.write('PSF CMAP CHEQ\n\n !NTITLE\n')
    ff.write("* \n")
    ff.write('%6d !NATOM\n' % natom);
    for i in range(ntriad):
        ii = i*4+1
        ff.write('%8d  %4s  %04d  %3s  O     0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % ( ii,segname,i,dname) );
        for j in range(3): # 3 arms of a triad
            ii = ii + 1
            ff.write('%8d  %4s  %04d  %3s  S     0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % (ii,segname,i,dname) );
    ff.write('\n%8d !NBOND\n' % nbond)
    cnt=0
    for i in range(ntriad):
        ii = i*4 +1
        for j in range(3):
            jj = ii + j + 1
            ff.write('%8d%8d' % (ii,jj))
            cnt = cnt + 1
            if cnt % 4 == 0 : ff.write('\n')
    return natom, nbond

def write_cm_psf(tva,tvb,cma,cmb,cha,chb):
    """ Write cm data as a psf """

    ff = open(fdir+"boc_cm.psf",'w')
    nframe=len(tva[0])
    natom = nframe * 6
    nbond = 0

    segname=['tva', 'tvb', 'cma', 'cmb', 'cha', ' chb']
    
    ff.write('PSF CMAP CHEQ\n\n !NTITLE\n')
    ff.write("* \n")
    ff.write('%6d !NATOM\n' % natom)

    for j in range(6):
        for i in range(nframe):
            jj= i+ j*nframe
            ff.write('%8d  %4s  %04d  %3s  CLA   0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % ( jj,segname[j],j,"BOC") );
    ff.write('\n%8d !NBOND\n' % nbond)
    ff.close()
    return

def write_cm_pdb(tva,tvb,cma,cmb,cha,chb):
    """ Write cm data as a pdb """    

    ff = open(fdir+"boc_cm.pdb",'w')
    nframe=len(tva[0])
    natom = nframe * 6

    segname=['tva', 'tvb', 'cma', 'cmb', 'cha', ' chb']
    ss=[tva[0],tvb[0],cma,cmb,cha,chb]
    
    occupancy = 0.5
    v = np.zeros(3); w = np.zeros(3)
    
    ff.write('MODEL\n');
    for j in range(6):
        for i in range(nframe):
            jj= i+ j*nframe
            v=ss[j][i]
            ff.write('ATOM  %5d  O   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                     (jj,"BOC",j,v[0],v[1],v[2],occupancy,segname[j]))
    ff.write('ENDMDL\n');
    ff.close()
    return

def write_pdb(segname,dname,triad):
    """ Write triad trajectory into ofname.pdb """
    ff = open(fdir+dname+".pdb",'w')
    ntriad =1 # use a single triad for now
    natom = ntriad*4
    occupancy = 0.5
    v = np.zeros(3); w = np.zeros(3)
    nframe = len(triad[0])
    for frame in range(nframe):
        ff.write('MODEL\n');
        for i in range(3):
            v[i] = triad[0][frame][i]
        ff.write('ATOM  %5d  O   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                 (1,dname,1,v[0],v[1],v[2],occupancy,segname))
        for j in range(3):
            if j == 0: tt = triad[1][frame]
            elif j == 1: tt = triad[2][frame]
            else: tt = triad[3][frame]
            jj = j + 2
            for k in range(3):
                w[k]=v[k]+15 * tt[k]
            ff.write('ATOM  %5d  S   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                     (jj,dname,1,w[0],w[1],w[2],occupancy,segname))
        ff.write('ENDMDL\n');

def write_boc_psf(ofname):
    """ Write boc psf """
    ff = open(fdir+'boc_'+ofname+".psf",'w')
    natom = 6
    nbond = 4
    ff.write('PSF CMAP CHEQ\n\n !NTITLE\n')
    ff.write("* \n")
    ff.write('%6d !NATOM\n' % natom);

    for i in range(2):
        ii = i*3+1
        ff.write('%8d  %4s  %04d  %3s  P     0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % ( ii,segi[i],1,vn[i]) );
        ii = i*3+2
        ff.write('%8d  %4s  %04d  %3s  P     0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % ( ii,segi[i],2,hn[i]) );
        ii = i*3+3
        ff.write('%8d  %4s  %04d  %3s  P     0   0.000000        1.0000               0   0.00000      0.000000E-02\n' % ( ii,segi[i],3,cn[i]) );

    ff.write('\n%8d !NBOND\n' % nbond)
    cnt=0
    for i in range(1,6):
        if i == 3: continue # no bond between atom 3 & 4
        ff.write('%8d%8d' % (i,i+1))
        #if i % 4 == 0 : ff.write('\n')
    return natom, nbond
    
def write_boc_pdb(molname,tva,tvb,cma,cmb,cha,chb):
    """ Write triad trajectory into ofname.pdb """
    ff = open(fdir+'boc_'+molname+".pdb",'w')
    occupancy = 0.5
    rv = np.zeros(3);  # var domain
    rc = np.zeros(3);  # const domain
    rh = np.zeros(3);  # hinge
    nframe = len(tva[0])
    for frame in range(nframe):
        ff.write('MODEL\n');
        for i in range(2):
            for j in range(3):
                if i == 0:
                    rv[j] = tva[0][frame][j];
                    rc[j] = cma[frame][j]; 
                    rh[j] = cha[frame][j];
                else: 
                    rv[j] = tvb[0][frame][j];
                    rc[j] = cmb[frame][j]; 
                    rh[j] = chb[frame][j];
            jj = i*3+1
            ff.write('ATOM  %5d  P   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                     (jj,vn[i],1,rv[0],rv[1],rv[2],occupancy,segi[i]))
            jj = i*3+2
            ff.write('ATOM  %5d  P   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                     (jj,hn[i],2,rh[0],rh[1],rh[2],occupancy,segi[i]))
            jj = i*3+3
            ff.write('ATOM  %5d  P   %3s  %04d    %8.3f%8.3f%8.3f%6.2f 0.00    %4s\n' % 
                     (jj,cn[i],3,rc[0],rc[1],rc[2],occupancy,segi[i]))
        ff.write('ENDMDL\n');

def write_boc_data(molname,tva,tvb,cma,cmb,cha,chb):
    """ Write triad trajectory into ofname.dat """
    ff = open(fdir+'boc_'+molname+".dat",'w')
    occupancy = 0.5
    rv = np.zeros(3);  # var domain
    rc = np.zeros(3);  # const domain
    rh = np.zeros(3);  # hinge
    nframe = len(tva[0])
    ff.write('# nframe= %5d\n' % (nframe))
    ff.write('# data order:\n')
    ff.write('# row 0-3: triad[0-3] row 4,5: r_hinge,r_const for triad_a\n')
    ff.write('# row 6-11: ditto for triad_b\n')
    for frame in range(nframe):
        ff.write('# %6d\n' %(frame))
        for j in range(4):
            ff.write('% 12.8e % 12.8e % 12.8e\n' \
                     % (tva[j][frame][0],tva[j][frame][1],tva[j][frame][2]))
        ff.write('% 12.8e % 12.8e % 12.8e\n' % (cha[frame][0],cha[frame][1],cha[frame][2]))
        ff.write('% 12.8e % 12.8e % 12.8e\n' % (cma[frame][0],cma[frame][1],cma[frame][2]))
        for j in range(4):
            ff.write('% 12.8e % 12.8e % 12.8e\n' \
                     % (tvb[j][frame][0],tvb[j][frame][1],tvb[j][frame][2]))
        ff.write('% 12.8e % 12.8e % 12.8e\n' % (chb[frame][0],chb[frame][1],chb[frame][2]))
        ff.write('% 12.8e % 12.8e % 12.8e\n' % (cmb[frame][0],cmb[frame][1],cmb[frame][2]))
    return

def write_cm_vmd(molname,tva,tvb,cma,cmb,cha,chb):
    """ Write cm data as a dot plot in vmd """
    ff = open(fdir+'cm_'+molname+".vmd",'w')

    q0=tva[0]
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    
    q0=tvb[0]
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    q0=cha
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    q0=chb
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    q0=cma
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    q0=cmb
    nframe = len(q0)
    ff.write('mol new\n')
    ff.write('graphics top color Cyan2\n')
    for frame in range(nframe):
        x=q0[frame][0]; y=q0[frame][1]; z=q0[frame][2]
        ff.write('graphics top sphere { %7.2f %7.2f %7.2f } radius .25\n'\
                 % (x,y,z))
    ff.write('\n')

    return

#######################################################

if __name__ == "__main__":

    tva=triad_gen(vn[0])
    tvb=triad_gen(vn[1])
    cma=read_cm(fdir+'cm_'+cn[0]+'.dat')
    cmb=read_cm(fdir+'cm_'+cn[1]+'.dat')
    cha=read_cm(fdir+'cm_'+hn[0]+'.dat')
    chb=read_cm(fdir+'cm_'+hn[1]+'.dat')

    write_psf(segi[0],vn[0])
    write_pdb(segi[0],vn[0],tva)
    write_psf(segi[1],vn[1])
    write_pdb(segi[1],vn[1],tvb)
    write_boc_psf(molname)
    write_boc_pdb(molname,tva,tvb,cma,cmb,cha,chb)
    write_boc_data(molname,tva,tvb,cma,cmb,cha,chb)
    write_cm_psf(tva,tvb,cma,cmb,cha,chb)
    write_cm_pdb(tva,tvb,cma,cmb,cha,chb)
    
    quit()
    
