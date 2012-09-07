import numpy as npy
import copy
import pylab
import pyn_fort_math as f


def average(a):

    leng=len(a)

    if leng > 0 :

        v,sigm=f.stats.average(a,leng)

    else :
        
        v,sigm=0.0,0.0

    return v,sigm



def histogram(a,bins=20,segment=None,delta_x=None,prec=None,norm=None,plot=True):

    leng=len(a)

    if norm==None:
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

    f.stats.histograma(opt_norm,opt_prec,opt_range,opt,a,bins,mmn,mmx,delta_x,prec,leng)

    h_x=copy.deepcopy(f.stats.histo_x)
    h_y=copy.deepcopy(f.stats.histo_y)
    f.stats.free_mem()
    if plot:
        pylab.plot(h_x,h_y,'ro-')

    return h_x,h_y


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
        h_x=copy.deepcopy(f.stats.histo_x)
        f.stats.free_mem()
        
        return h_x

    else:

        tray_bins=f.stats.binning(opt_prec,opt_range,opt,traj,bins,mmn,mmx,delta_x,prec,len(traj))

        h_x=copy.deepcopy(f.stats.histo_x)
    
        f.stats.free_mem()

        return h_x,tray_bins


