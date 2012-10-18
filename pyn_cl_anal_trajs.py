import pyn_fort_anal_trajs as ftrajs
from pyn_cl_net import *
import pyn_math as pyn_math
from numpy import *
import copy as ccopy

## traj: frame,num_part,dim

class kinetic_1D_analysis():

    def __init__ (self,traject=None,column=None,num_particles=None,dimensions=None,verbose=True):

        self.file_name=''
        self.file_column=None
        self.traj=None
        self.dimensions=0
        self.num_particles=0
        self.num_frames=0
        self.traj_nodes=None
        self.traj_clusters=None
        self.network_nodes=None
        self.network_clusters=None

        if traject==None:
            print 'Trajectory needed (name of file or array)'
            return
        
        if type(column) in [tuple,list]:
            self.num_particles=len(column)
            self.traj=[[] for ii in range(len(self.num_particles))]
        elif type(column) in [int]:
            self.num_particles=1
            self.traj=[]
            self.traj_nodes=[]
            self.traj_clusters=[]
        elif column in ['ALL','All','all']:
            fff=open(traject,'r')
            line=fff.readline().split()
            self.num_particles=len(line)
            fff.close()
            self.traj=[[] for ii in range(len(self.num_particles))]
            self.traj_nodes=[[] for ii in range(len(self.num_particles))]
            self.traj_clusters=[[] for ii in range(len(self.num_particles))]

        if type(traject) in [str]:

            if type(column) in [tuple,list]:
                self.num_particles=len(column)
                self.traj=[[] for ii in range(len(self.num_particles))]
            elif type(column) in [int]:
                self.num_particles=1
                self.traj=[]
                self.traj_nodes=[]
                self.traj_clusters=[]
            elif column in ['ALL','All','all']:
                fff=open(traject,'r')
                line=fff.readline().split()
                self.num_particles=len(line)
                fff.close()
                self.traj=[[] for ii in range(len(self.num_particles))]
                self.traj_nodes=[[] for ii in range(len(self.num_particles))]
                self.traj_clusters=[[] for ii in range(len(self.num_particles))]

            self.file_name=traject
            self.file_column=column
            
            fff=open(traject,'r')
            for line in fff:
                line=line.split()
                self.traj.append(float(line[self.file_column]))
 
            fff.close()

        elif type(traject) in [list,tuple,ndarray]:
            self.traj=traject
            if num_particles==None and dimensions==None:
                self.num_particles=1
                self.dimensions=1
                self.traj_nodes=[]
                self.traj_clusters=[]

        if type(self.traj) not in [ndarray]:
            self.traj=array(self.traj,order="Fortran")

        self.num_frames=self.traj.shape[-1]
        self.dimensions=1

        if verbose:
            if self.num_particles==1:
                print '# Trajectory loaded:',self.num_frames,'time steps'
            else:
                print '# Trajectory loaded:',self.num_particles,'particles,',self.num_frames,'time steps'

    def histogram(self,bins=20,segment=None,delta_x=None,norm=False,cumul=False):

        if self.num_particles==1:
            return pyn_math.histogram(self.traj,bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)
        else:
            return pyn_math.histogram(self.traj[0],bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)

    def life_time(self,traj=None,state=None,segment=None,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1

        if traj == None:

            if type(state) in [int,float]:
                num_states=1
                state=[state]
            elif type(state) in [list,tuple]:
                num_states=len(state)

            traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
            lt_mean=ftrajs.aux.life_time_dist(opt_norm,traj_inp,state,self.num_frames,self.num_particles,self.dimensions,num_states)


        elif traj in ['CLUSTERS','Clusters','clusters']:

            if type(state) in [int]:
                num_states=1
                state=[state]
            elif type(state) in [list,tuple]:
                num_states=len(state)

            traj_inp=standard_traj(self.traj_clusters,num_parts=self.num_particles,dims=self.dimensions)
            lt_mean=ftrajs.aux.life_time_dist(opt_norm,traj_inp,state,self.num_frames,self.num_particles,self.dimensions,num_states)

        elif traj in ['NODES','Nodes','nodes']:

            if type(state) in [int]:
                num_states=1
                state=[state]
            elif type(state) in [list,tuple]:
                num_states=len(state)

            traj_inp=standard_traj(self.traj_nodes,num_parts=self.num_particles,dims=self.dimensions)
            lt_mean=ftrajs.aux.life_time_dist(opt_norm,traj_inp,state,self.num_frames,self.num_particles,self.dimensions,num_states)

    
        lt_dist=ccopy.deepcopy(ftrajs.aux.distrib)
        lt_x=ccopy.deepcopy(ftrajs.aux.distrib_x)
        ftrajs.aux.free_distrib()

        if verbose:
            print '# Mean life time:', lt_mean,'frames.'

        if mean:
            return lt_mean
        else:
            return lt_x, lt_dist

    def first_passage_time (self,traj=None,from_state=None,from_segment=None,to_state=None,to_segment=None,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_to_segment=0
        opt_from_segment=0
        opt_to_state=0
        opt_from_state=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1

        if from_state!=None:   
            opt_from_state=1
        else:
            from_state=0

        if from_segment!=None: 
            opt_from_segment=1
        else:
            from_segment=[0.0,0.0]

        if to_state!=None:     
            opt_to_state=1
        else:
            to_state=0

        if to_segment!=None:   
            opt_to_segment=1
        else:
            to_segment=[0.0,0.0]
        
        if opt_to_segment==0 and opt_to_state==0:
            print '# the input variable to_state or to_segment is needed'
            pass

        if traj == None:
            traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=standard_traj(self.traj_clusters,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=standard_traj(self.traj_nodes,num_parts=self.num_particles,dims=self.dimensions)
        else:
            print '# A readable traj is needed'
            pass

        if type(from_state) in [int,float]:
            from_num_states=1
            from_state=[from_state]
        elif type(state) in [list,tuple]:
            from_num_states=len(from_state)

        if type(to_state) in [int,float]:
            to_num_states=1
            to_state=[to_state]
        elif type(to_state) in [list,tuple]:
            to_num_states=len(to_state)


        fpt_mean=ftrajs.aux.fpt_dist(opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,self.num_frames,self.num_particles,self.dimensions,\
                                                        from_num_states,to_num_states)

        fpt_dist=ccopy.deepcopy(ftrajs.aux.distrib)
        fpt_x=ccopy.deepcopy(ftrajs.aux.distrib_x)
        ftrajs.aux.free_distrib()

        if verbose:
            print '# Mean first passage time:', fpt_mean,'frames.'

        if mean:
            return fpt_mean
        else:
            return fpt_x, fpt_dist


    def first_committed_passage_time (self,traj=None,states=None,segments=None,commitment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_segments=0
        opt_states=0
        opt_noreturn=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_noreturn=1

        if states!=None:   
            opt_states=1
            segments=[[0,0]]
            num_segments=1
        else:
            opt_segments=0
            states=[0]
            num_states=1

        if opt_segments==0 and opt_states==0:
            print '# the input variable states or segments is needed'
            pass

        if traj == None:
            traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=standard_traj(self.traj_clusters,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=standard_traj(self.traj_nodes,num_parts=self.num_particles,dims=self.dimensions)
        else:
            print '# A readable traj is needed'
            pass

        if type(states) in [int,float]:
            num_states=1
            states=[states]
        elif type(states) in [list,tuple]:
            num_states=len(states)

        if opt_segments:
            num_segments=len(segments)

        if commitment==None:
            if opt_segments:
                num_commits=num_segments
            else:
                num_commits=num_states
            commitment=[True for ii in range(num_commits)]
        else:
            num_commits=len(commitment)

        commitment_in=[int(ii) for ii in commitment]
        fcpt_mean=ftrajs.aux.fcpt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,self.num_frames,self.num_particles,self.dimensions,\
                                                        num_states,num_segments,num_commits)

        fcpt_dist=ccopy.deepcopy(ftrajs.aux.distrib)
        fcpt_x=ccopy.deepcopy(ftrajs.aux.distrib_x)
        ftrajs.aux.free_distrib()

        if verbose:
            print '# Mean first passage time:', fcpt_mean,'frames.'

        if mean:
            return fcpt_mean
        else:
            return fcpt_x, fcpt_dist


    def trip_time (self,traj=None,from_state=None,from_segment=None,to_state=None,to_segment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_no_return=0
        opt_to_segment=0
        opt_from_segment=0
        opt_to_state=0
        opt_from_state=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_no_return=1

        if from_state!=None:   
            opt_from_state=1
        else:
            from_state=0

        if from_segment!=None: 
            opt_from_segment=1
        else:
            from_segment=[0.0,0.0]

        if to_state!=None:     
            opt_to_state=1
        else:
            to_state=0

        if to_segment!=None:   
            opt_to_segment=1
        else:
            to_segment=[0.0,0.0]
        
        if opt_to_segment==0 and opt_to_state==0:
            print '# the input variable to_state or to_segment is needed'
            pass

        if traj == None:
            traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=standard_traj(self.traj_clusters,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=standard_traj(self.traj_nodes,num_parts=self.num_particles,dims=self.dimensions)
        else:
            print '# A readable traj is needed'
            pass

        if type(from_state) in [int,float]:
            from_num_states=1
            from_state=[from_state]
        elif type(state) in [list,tuple]:
            from_num_states=len(from_state)

        if type(to_state) in [int,float]:
            to_num_states=1
            to_state=[to_state]
        elif type(to_state) in [list,tuple]:
            to_num_states=len(to_state)


        tt_mean=ftrajs.aux.tt_dist(opt_norm,opt_no_return,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, \
                                                        from_state,from_segment,to_state,to_segment, \
                                                        traj_inp,self.num_frames,self.num_particles,self.dimensions,\
                                                        from_num_states,to_num_states)

        tt_dist=ccopy.deepcopy(ftrajs.aux.distrib)
        tt_x=ccopy.deepcopy(ftrajs.aux.distrib_x)
        ftrajs.aux.free_distrib()

        if verbose:
            print '# Mean first passage time:', tt_mean,'frames.'

        if mean:
            return tt_mean
        else:
            return tt_x, tt_dist

    def committed_trip_time (self,traj=None,states=None,segments=None,commitment=None,no_return=False,mean=False,norm=False,verbose=False):

        opt_mean=0
        opt_norm=0
        opt_segments=0
        opt_states=0
        opt_noreturn=0

        if (mean):
            opt_mean=1
        if (norm):
            opt_norm=1
        if (no_return):
            opt_noreturn=1

        if states!=None:   
            opt_states=1
            segments=[[0,0]]
            num_segments=1
        else:
            opt_segments=0
            states=[0]
            num_states=1

        if opt_segments==0 and opt_states==0:
            print '# the input variable states or segments is needed'
            pass

        if traj == None:
            traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['CLUSTERS','Clusters','clusters']:
            traj_inp=standard_traj(self.traj_clusters,num_parts=self.num_particles,dims=self.dimensions)
        elif traj in ['NODES','Nodes','nodes']:
            traj_inp=standard_traj(self.traj_nodes,num_parts=self.num_particles,dims=self.dimensions)
        else:
            print '# A readable traj is needed'
            pass

        if type(states) in [int,float]:
            num_states=1
            states=[states]
        elif type(states) in [list,tuple]:
            num_states=len(states)

        if opt_segments:
            num_segments=len(segments)

        if commitment==None:
            if opt_segments:
                num_commits=num_segments
            else:
                num_commits=num_states
            commitment=[True for ii in range(num_commits)]
        else:
            num_commits=len(commitment)

        commitment_in=[int(ii) for ii in commitment]
        ctt_mean=ftrajs.aux.ctt_dist(opt_norm,opt_noreturn,opt_states,opt_segments, states,segments,commitment_in,\
                                                        traj_inp,self.num_frames,self.num_particles,self.dimensions,\
                                                        num_states,num_segments,num_commits)

        ctt_dist=ccopy.deepcopy(ftrajs.aux.distrib)
        ctt_x=ccopy.deepcopy(ftrajs.aux.distrib_x)
        ftrajs.aux.free_distrib()

        if verbose:
            print '# Mean first passage time:', ctt_mean,'frames.'

        if mean:
            return ctt_mean
        else:
            return ctt_x, ctt_dist


    def kinetic_network(self,traj=None,verbose=False):

        if traj in ['CLUSTERS','Clusters','clusters']:
            if type(self.traj_clusters) not in [ndarray]:
                self.traj_clusters=array(self.traj_clusters,order="Fortran")
            self.network_clusters=kinetic_network(self.traj_clusters,ranges=[self.traj_clusters.min(),self.traj_clusters.max()],verbose=verbose)
            pass

        elif traj in ['NODES','Nodes','nodes']:
            if type(self.traj_nodes) not in [ndarray]:
                self.traj_nodes=array(self.traj_nodes,order="Fortran")
            self.network_nodes=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=verbose)
            pass

        else:
            print 'traj= clusters or nodes'
            pass



    def berezovska2012(self,window=None,ksi=0.5,granularity=1.2,bins=20,segment=None,delta_x=None,clusters=True,verbose=False):

        if segment==None:
            opt_range=0
            mmx=self.traj.max()
            mmn=self.traj.min()
        else:
            opt_range=1
            mmn=segment[0]
            mmx=segment[1]

        if delta_x!=None:
            opt=1
        else:
            delta_x=1.0 # Its given by gannas function
            opt=2


        traj_inp=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
        self.traj_nodes=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,traj_inp,ksi,window,self.num_particles,self.num_frames)

        self.network=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=False)

        if clusters:

            self.network.symmetrize(new=False,verbose=verbose)
            self.network.mcl(granularity=granularity,pruning=True,verbose=verbose)

            if self.num_particles==1:
                self.traj_clusters=[]
                for ii in self.traj_nodes:
                    self.traj_clusters.append(self.network.node[ii].cluster)
            else:
                self.traj_cluster=[[] for ii in range(len(self.num_particles))]
                for jj in self.num_particles:
                    for ii in self.traj_nodes[jj,:]:
                        self.traj_clusters[jj].append(self.network.node[ii].cluster)
             
            self.traj_clusters=array(self.traj_clusters,order="Fortran")



    def fpt_rao(self,traj=None,target=0,norm=False):
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

    def rao(self,traj=None,window=None,separators=None):

        if traj==None or window==None or separators==None:
            print 'Not enough input variables.'
            pass
         
        traj_i=traj
        if type(traj_i) not in [ndarray]:
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




#def standard_traj(traj,num_parts,dims):
# 
#    if type(traj) not in [numpy.ndarray]:
#        traj=array(traj,order='Fortran')
# 
#    if len(traj.shape)==3:
#        return traj
#    elif len(traj.shape)==1:
#        return [[traj]]
#    else:
#        if num_parts==1:
#            return [traj]
#        else:  #dims==1
            


