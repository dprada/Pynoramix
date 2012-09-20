import pyn_fort_anal_trajs as ftrajs
from numpy import *

def ganna (traj=None,window=None,ksi=0.5,bins=20,segment=None,delta_x=None):    

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

    traj_i=traj
    if type(traj_i) not in ['numpy.ndarray']:
        traj_i=array(traj_i,order='Fortran')

    if len(traj_i.shape)==1:
        traj_i=array([traj_i],order='Fortran')

    num_parts=traj_i.shape[0]
    frames=traj_i.shape[1]
    clust_trajs=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,traj_i,ksi,window,num_parts,frames)

    if num_parts==1:
        return clust_trajs[0]
    else:
        return clust_trajs

def ganna2 (traj=None,window=None,ksi=0.5,bins=20,segment=None,delta_x=None):    

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

    traj_i=traj
    if type(traj_i) not in ['numpy.ndarray']:
        traj_i=array(traj_i,order='Fortran')

    if len(traj_i.shape)==1:
        traj_i=array([traj_i],order='Fortran')

    num_parts=traj_i.shape[0]
    frames=traj_i.shape[1]
    clust_trajs=ftrajs.aux.ganna2(opt_range,opt,bins,mmn,mmx,delta_x,traj_i,ksi,window,num_parts,frames)

    if num_parts==1:
        return clust_trajs[0]
    else:
        return clust_trajs


def fpt_rao(traj=None,target=0,norm=False):
    distrib={}
    ant=traj[-1]
    counter=0
    for ii in range(1,len(traj)):
        act=traj[-(ii+1)]
        if ant==target:
            counter=0
        counter+=1
        if act!=target:
            try:
                distrib[counter]+=1
            except:
                distrib[counter]=1
        ant=act
    xx=distrib.keys()
    yy=distrib.values()
    if norm:
        lll=sum(yy[:])
        for ii in range(len(yy)):
            yy[ii]=(yy[ii]*1.0)/(lll*1.0)
    return xx,yy

def rao(traj=None,window=None,separators=None):

    if traj==None or window==None or separators==None:
        print 'Not enough input variables.'
        pass

    traj_i=traj
    if type(traj_i) not in ['numpy.ndarray']:
        traj_i=array(traj_i,order='Fortran')

    if len(traj_i.shape)==1:
        traj_i=array([traj_i],order='Fortran')

    num_parts=traj_i.shape[0]
    frames=traj_i.shape[1]
    salida=ftrajs.aux.rao_stat_1(window,traj,separators,num_parts,frames,len(separators))

    if num_parts==1:
        return salida[0]
    else:
        return salida
