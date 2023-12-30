"""
anal_triad.py: read triad data and make various measurements

  tv[ab]: centroid and triad for V domains
  cm[ab]: center of mass of C domains
  ch[ab]: center of mass of hinges
"""
import sys
import numpy as np

from name_def import get_names
segi, vn,cn,hn,molname, fdir = get_names()

#######################################################
def read_boc_data(molname):
    """ Read triad trajectory """
    ff = open(fdir+'boc_'+molname+".dat",'r')
    nframe=-1
    ifrm=-1 # number of frames read
    tva=[[] for i in range(4) ]; tvb=[[] for i in range(4)]
    cma=[]; cmb=[]; cha=[]; chb=[]
    while 1:
        line=ff.readline()
        if not line: break
        a=line.split()
        if a[0]=='#':
            if a[1]=='nframe=': nframe=int(a[2])
            if nframe>0:
                if len(a)==2: # frame number read
                    ifrm=int(a[1])
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tva[0].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tva[1].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tva[2].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tva[3].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    cha.append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    cma.append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tvb[0].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tvb[1].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tvb[2].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    tvb[3].append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    chb.append(a.tolist())
                    line1=ff.readline()
                    a=np.asfarray(np.array(line1.split()),float)
                    cmb.append(a.tolist())
    print('nframe= ', nframe)
    return nframe,tva,tvb,cha,chb,cma,cmb

def angle_hinge(nframe,tva,tvb,cha,chb,cma,cmb):
    """ Measure hinge angles for a and b chains """
    phi=np.zeros((2,nframe))
    for i in range(nframe):
        da1=np.subtract(cma[i],cha[i])
        db1=np.subtract(cmb[i],chb[i]);
        da2=np.subtract(tva[0][i],cha[i])
        db2=np.subtract(tvb[0][i],chb[i])
        pa=np.dot(da1,da2); pb=np.dot(db1,db2)
        na1=np.linalg.norm(da1)
        na2=np.linalg.norm(da2)
        nb1=np.linalg.norm(db1)
        nb2=np.linalg.norm(db2)
        pa=pa/(na1*na2); pb=pb/(nb1*nb2)
        phi[0][i]=np.arccos(pa)
        phi[1][i]=np.arccos(pb)
        
    return phi

def angle_twist(nframe,tva,tvb,cha,chb,cma,cmb):
    """ Measure twist angle between V domains """
    psi=np.zeros(nframe)
    for i in range(nframe):
        de0=np.subtract(tva[0][i],tvb[0][i]) # centroid-centroid vector betwen V domains
        de0=normalize(de0)
        e3a=tva[3][i]; e3a = e3a - np.dot(e3a,de0) * de0 ; e3a=normalize(e3a)
        e3b=tvb[3][i]; e3b = e3b - np.dot(e3b,de0) * de0 ; e3b=normalize(e3b)
        rdum=np.dot(e3a,e3b)
        psi[i]=np.arccos(rdum)

    return psi

def angle_dihe(nframe,a,b,c,d):
    """ Measure dihedral angle for 4 atoms a,b,c,d.
        a-d: cm position vectors.
        sign: + if \vec{cd} rotates counterclockwise relative to \vec{ab} about \vec{bc}
    """
    psi=np.zeros(nframe)
    for i in range(nframe):
        ba=np.subtract(a[i],b[i]);        ba=normalize(ba)
        bc=np.subtract(c[i],b[i]);        bc=normalize(bc)
        cd=np.subtract(d[i],c[i]);        cd=normalize(cd)
        
        ea = ba - np.dot(ba,bc) * bc ; ea=normalize(ea)
        ed = cd - np.dot(cd,bc) * bc ; ed=normalize(ed)
        rdum=np.dot(ea,ed)
        abcd=np.cross(ba,cd)
        if np.dot(abcd,bc)>0.: rdum1=1.
        else: rdum1=-1.
        psi[i]=rdum1*np.arccos(rdum)        

    return psi

def dist_boc(nframe,tva,tvb,cha,chb,cma,cmb):
    """ Measure distances """
    dist=np.zeros((2,nframe))
    for i in range(nframe):
        dc=np.subtract(cma[i],cmb[i])
        dv=np.subtract(tva[0][i],tvb[0][i])
        dist[0][i]=np.linalg.norm(dv)
        dist[1][i]=np.linalg.norm(dc)
        
    return dist

def normalize(a):
    """ normalize a vector """
    return a/ np.linalg.norm(a)

def write_data(nframe,phi,dist,psi,psi_a,psi_b):
    """ Write data """
    ff = open(fdir+'meas_'+molname+".dat",'w')
    ff.write('# phi_a(deg) phi_b(deg)  d(Va-Vb)  d(Ca-Cb) psi_VV(deg)   psi_a   psi_b\n')
    for i in range(nframe):
        a0=phi[0][i]*180./np.pi
        a1=phi[1][i]*180./np.pi
        a2=psi[i]*180./np.pi
        a3=psi_a[i]*180./np.pi
        a4=psi_b[i]*180./np.pi
        ff.write('% 10.5f % 10.5f % 10.5f % 10.5f %10.5f %10.5f %10.5f\n' % \
                 (a0,a1,dist[0][i],dist[1][i],a2,a3,a4))
    return

#######################################################

if __name__ == "__main__":

    nframe,tva,tvb,cha,chb,cma,cmb = read_boc_data(molname)
    phi=angle_hinge(nframe,tva,tvb,cha,chb,cma,cmb)
    psi=angle_twist(nframe,tva,tvb,cha,chb,cma,cmb)
    dist=dist_boc(nframe,tva,tvb,cha,chb,cma,cmb)

    # center of mass of Cab
    cm_avg=0.5*(np.asarray(cma)+np.asarray(cmb))
    # dihedral angles
    psi_a=angle_dihe(nframe,chb,cm_avg,cha,tva[0]) #  Cb-Cab-ha-Va
    psi_b=angle_dihe(nframe,cha,cm_avg,chb,tvb[0]) #  Ca-Cab-hb-Vb

    write_data(nframe,phi,dist,psi,psi_a,psi_b)
    quit()
    
