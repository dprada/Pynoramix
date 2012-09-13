import pyn_fort_anal_trajs as ftrajs

def ganna (traj=None,window=None,ksi=0.5,bins=20,segment=None,delta_x=None,prec=None):    
    leng=len(traj)
    opt_norm=1 #(norm=True)
    opt_cumul=1 #(cumul=True)
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
    clust_trajs=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,traj,ksi,window,leng)
    return clust_trajs


#ftrajs.aux.rao_stat_1(tw,result[:,ii,0],bond,len(result),len(bond))
