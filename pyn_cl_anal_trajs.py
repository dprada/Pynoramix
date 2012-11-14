import pyn_fort_anal_trajs as ftrajs
from pyn_cl_net import *
import pyn_math as pyn_math
from numpy import *
import copy as ccopy

## traj: frame,num_part,dim

class kinetic_1D_analysis():

    def __init__ (self,traject=None,columns=None,particles=None,dimensions=None,verbose=True):

        self.file_name=''
        self.file_columns=None
        self.traj=None
        self.dimensions=dimensions
        self.particles=particles
        self.frames=0
        self.traj_nodes=None
        self.traj_clusters=None
        self.network_nodes=None
        self.network_clusters=None
        self.__type_nodes__=None
        self.__type_clusters__=None
        self.__offset__=None

        if traject==None:
            print 'Trajectory needed (name of file or array)'
            return
        
        if type(columns) in [int]:
            self.particles=1
            self.dimensions=1
            columns=[columns]
        elif type(columns) in [tuple,list]:
            if (len(columns))==1:
                self.particles=1
                self.dimensions=1
            else:
                if particles==None and dimensions==None:
                    print '# Please, make use of the variables "particles" or/and "dimensions":'
                    print '#   traj:=[100 frames, 3 dimensions] --> "particles=1" or/and "dimensions=3"'
                    print '#   traj:=[100 frames, 8  particles] --> "particles=8" or/and "dimensions=1"'
                    return
                elif particles==1 or dimensions>=1:
                    self.particles=1
                elif particles>1 or dimensions==1:
                    self.dimensions=1

        elif columns in ['ALL','All','all']:
            fff=open(traject,'r')
            line=fff.readline().split()
            nn=len(line)
            fff.close()
            columns=[ii for ii in range(nn)]
            if nn>1:
                if particles==None and dimensions==None:
                    print '# Please, make use of the variables "particles" or/and "dimensions":'
                    print '#   traj:=[100 frames, 3 dimensions] --> "particles=1" or/and "dimensions=3"'
                    print '#   traj:=[100 frames, 8  particles] --> "particles=8" or/and "dimensions=1"'
                    return
                elif particles==1 or dimensions>=1:
                    self.particles=1
                    self.dimensions=nn
                elif particles>1 or dimensions==1:
                    self.dimensions=1
                    self.particles=nn
            else:
                self.dimensions=1
                self.particles=1

        if type(traject) in [str]:

            self.file_name=traject
            self.file_column=columns
            
            fff=open(traject,'r')
            for line in fff:
                line=line.split()
                self.traj.append([float(line[ii]) for ii in columns])
 
            fff.close()

        if type(traject) in [list,tuple,ndarray]:
            
            self.traj=standard_traj(traject,particles=self.particles,dimensions=self.dimensions)

        self.frames=self.traj.shape[0]
        self.particles=self.traj.shape[1]
        self.dimensions=self.traj.shape[2]

        if verbose:
            print '# Loaded:',self.frames,'frames,',self.particles,'particles,',self.dimensions,'dimensions.'


    def histogram1D(self,dimension=None,node=None,cluster=None,bins=20,segment=None,delta_x=None,norm=False,cumul=False):

        if cluster==None and node==None:
            return pyn_math.histogram1D(self.traj,bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)

        if cluster!=None:

            if self.num_particles==1:
                
                if self.__type_clusters__=='berezovska2012':

                    traj_inp=[]
                    tw=(self.num_frames-len(self.traj_clusters))/2
                    for ii in range(len(self.traj_clusters)):
                        if (self.traj_clusters[ii]==cluster):
                            traj_inp.append(self.traj[ii+tw])

                return pyn_math.histogram(traj_inp,bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)

            else:
                print 'Not included yet'
                pass

        if node!=None:

            print 'Not included yet'
            pass

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
                                                        traj_inp,traj_inp.shape[0],traj_inp.shape[1],traj_inp.shape[2],\
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


    def kinetic_network(self,traj=None,ranges=None,verbose=False):

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
            self.traj=standard_traj(self.traj,num_parts=self.num_particles,dims=self.dimensions)
            if ranges==None:
                ranges=build_ranges(self.traj)
            else:
                ranges=standard_ranges(ranges)
            self.network_nodes,self.traj_nodes=kinetic_network(self.traj,ranges=ranges,traj_out=True,verbose=verbose)                
            pass



    def berezovska2012(self,window=None,ksi=0.5,granularity=1.2,bins=20,segment=None,delta_x=None,clusters=True,verbose=False):

        ref_max=self.traj.max()
        ref_min=self.traj.min()
        rv_min=0
        rv_max=0

        if segment==None:
            opt_range=0
            mmx=ref_max
            mmn=ref_min
        else:
            opt_range=1
            mmn=segment[0]
            mmx=segment[1]
            if mmn>ref_min:
                rv_min=1
            if mmx<ref_max:
                rv_max=1

        if delta_x!=None:
            opt=1
        else:
            delta_x=1.0 # Its given by gannas function
            opt=2

        if self.dimensions!=1:
            print '# Method not implemented yet for more than 1D.'
            return

        if verbose:
            if rv_min:
                print '# Extra node for x <', mmn
            if rv_max:
                print '# Extra node for x >', mmx

        self.traj_nodes=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,rv_min,rv_max,self.traj,ksi,window,self.num_particles,self.num_frames)

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

        self.__type_clusters__='berezovska2012'




#    def rao(self,traj=None,window=None,separators=None):
# 
#        if traj==None or window==None or separators==None:
#            print 'Not enough input variables.'
#            pass
#         
#        traj_i=traj
#        if type(traj_i) not in [ndarray]:
#            traj_i=array(traj_i,order='Fortran')
#         
#        if len(traj_i.shape)==1:
#            traj_i=array([traj_i],order='Fortran')
#         
#        num_parts=traj_i.shape[0]
#        frames=traj_i.shape[1]
#        salida=ftrajs.aux.rao_stat_1(window,traj,separators,num_parts,frames,len(separators))
#         
#        if num_parts==1:
#            return salida[0]
#        else:
#            return salida




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
            


