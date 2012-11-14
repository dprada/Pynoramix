import numpy as npy
import copy as ccopy
import pyn_fort_math as fstat
try:
    import pylab
    wpylab=1
except:
    wpylab=0

def average(a):

    leng=len(a)

    if leng > 0 :

        v,sigm=fstat.stats.average(a,leng)

    else :
        
        v,sigm=0.0,0.0

    return v,sigm



def histogram1D(traj,bins=20,segment=None,delta_x=None,norm=False,cumul=False,plot=False):

    if type(traj) not in [npy.ndarray]:
        traj=npy.array(traj,dtype=float,order='F')

    if len(traj.shape)==1:
        traj.resize(1,traj.shape[0])

    part=traj.shape[0]
    leng=traj.shape[1]

    if not npy.isfortran(traj):
        traj=npy.array(traj,dtype=float,order='F')

    if norm==False:
        opt_norm=0
    else:
        opt_norm=1

    if segment==None:
        opt_range=0
        mmx=traj.max()
        mmn=traj.min()
    else:
        opt_range=1
        mmn=segment[0]
        mmx=segment[1]

    if delta_x!=None:
        opt=1
    else:
        delta_x=1.0
        opt=2

    if cumul:
        opt_cumul=1
    else:
        opt_cumul=0

    fstat.stats.histogram1d(opt_norm,opt_range,opt,opt_cumul,traj,bins,mmn,mmx,delta_x,part,leng)

    h_x=ccopy.deepcopy(fstat.stats.histo_x)
    h_y=ccopy.deepcopy(fstat.stats.histo_y)
    fstat.stats.free_mem()
    if plot and wpylab:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y

def histogram2D(traj,bins=[20,20],segment=None,delta_x=None,prec=None,norm=False,plot=False):

    leng=len(traj)

    if norm==False:
        opt_norm=0
    else:
        opt_norm=1

    if prec==None:
        opt_prec=0
        prec=1.0
    else:
        opt_prec=1

    if segment==None:
        opt_range=0
        mmx0=max(traj[:,0])
        mmx1=max(traj[:,1])
        mmn0=min(traj[:,0])
        mmn1=min(traj[:,1])
    else:
        opt_range=1
        mmn0=segment[0][0]
        mmn1=segment[1][0]
        mmx0=segment[0][1]
        mmx1=segment[1][1]

    if delta_x!=None:
        opt=1
    else:
        delta_x=[1.0,1.0]
        opt=2

    fstat.stats.histograma_2d(opt_norm,opt_prec,opt_range,opt,traj,bins,[mmn0,mmn1],[mmx0,mmx1],delta_x,prec,leng)

    h_x=ccopy.deepcopy(fstat.stats.histo_x)
    h_y=ccopy.deepcopy(fstat.stats.histo_y)
    h_z=ccopy.deepcopy(fstat.stats.histo_z)
    fstat.stats.free_mem()
    if plot and wpylab:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y,h_z

def binning(traj=None,bins=20,segment=None,delta_x=None,prec=None):


    if segment==None:
        opt_range=0
        mmx=0.0
        mmn=0.0
    else:
        opt_range=1
        mmn=segment[0]
        mmx=segment[1]

    if delta_x==None:
        opt_delta_x=0
        delta_x=1.0
    else:
        opt_delta_x=1

    if traj==None:

        o_delta_x=fstat.stats.binning_x(opt_range,opt,bins,mmn,mmx,delta_x)
        h_x=ccopy.deepcopy(fstat.stats.histo_x)
        fstat.stats.free_mem()
        
        return h_x

    else:

        tray_bins=fstat.stats.binning(opt_range,opt_delta_x,traj,bins,mmn,mmx,delta_x,len(traj))

        h_x=ccopy.deepcopy(fstat.stats.histo_x)
    
        fstat.stats.free_mem()

        return h_x,tray_bins


