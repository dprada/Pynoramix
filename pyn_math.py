import numpy as npy
import copy as ccopy
import pylab
import pyn_fort_math as f


def average(a):

    leng=len(a)

    if leng > 0 :

        v,sigm=f.stats.average(a,leng)

    else :
        
        v,sigm=0.0,0.0

    return v,sigm



def histogram(traj,bins=20,segment=None,delta_x=None,prec=None,norm=False,cumul=False,plot=False):

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
        mmx=max(traj)
        mmn=min(traj)
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

    f.stats.histograma(opt_norm,opt_prec,opt_range,opt,opt_cumul,traj,bins,mmn,mmx,delta_x,prec,leng)

    h_x=ccopy.deepcopy(f.stats.histo_x)
    h_y=ccopy.deepcopy(f.stats.histo_y)
    f.stats.free_mem()
    if plot:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y

def histogram2d(traj,bins=[20,20],segment=None,delta_x=None,prec=None,norm=False,plot=False):

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

    f.stats.histograma_2d(opt_norm,opt_prec,opt_range,opt,traj,bins,[mmn0,mmn1],[mmx0,mmx1],delta_x,prec,leng)

    h_x=ccopy.deepcopy(f.stats.histo_x)
    h_y=ccopy.deepcopy(f.stats.histo_y)
    h_z=ccopy.deepcopy(f.stats.histo_z)
    f.stats.free_mem()
    if plot:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y,h_z

def binning(traj=None,bins=20,segment=None,delta_x=None,prec=None):


    if prec==None:
        opt_prec=0
        prec=1.0
    else:
        opt_prec=1

    if segment==None:
        opt_range=0
        mmx=max(a)
        mmn=min(a)
    else:
        opt_range=1
        mmn=segment[0]
        mmx=segment[1]

    if delta_x!=None:
        opt=1
    else:
        delta_x=1.0
        opt=2

    if traj==None:

        o_delta_x=f.stats.binning_x(opt_range,opt,bins,mmn,mmx,delta_x,prec)
        h_x=ccopy.deepcopy(f.stats.histo_x)
        f.stats.free_mem()
        
        return h_x

    else:

        tray_bins=f.stats.binning(opt_prec,opt_range,opt,traj,bins,mmn,mmx,delta_x,prec,len(traj))

        h_x=ccopy.deepcopy(f.stats.histo_x)
    
        f.stats.free_mem()

        return h_x,tray_bins


