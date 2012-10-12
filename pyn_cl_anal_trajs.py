import pyn_fort_anal_trajs as ftrajs
from pyn_cl_net import *
import pyn_math as pyn_math
from numpy import *

class kinetic_1D_analysis():

    def __init__ (self,traject=None,column=None,verbose=True):

        self.file_name=''
        self.file_column=None
        self.traj=None
        self.dimensions=1
        self.num_particles=0
        self.length=0
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
            self.file_name=traject
            self.file_column=column
            
            fff=open(traject,'r')
            for line in fff:
                line=line.split()
                self.traj.append(float(line[self.file_column]))
 
            fff.close()

        if type(self.traj) not in ['numpy.ndarray']:
            self.traj=array(self.traj,order="Fortran")

        self.length=self.traj.shape[-1]

        if verbose:
            if self.num_particles==1:
                print '# Trajectory loaded:',self.length,'time steps'
            else:
                print '# Trajectory loaded:',self.num_particles,'particles,',self.length,'time steps'

    def histogram(self,bins=20,segment=None,delta_x=None,norm=False,cumul=False):

        if self.num_particles==1:
            return pyn_math.histogram(self.traj,bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)
        else:
            return pyn_math.histogram(self.traj[0],bins=bins,segment=segment,delta_x=delta_x,norm=norm,cumul=cumul,plot=False)

    def life_time(self,traj=None,state=None,mean=False,verbose=False):

        if traj == None:
            pass

        elif traj in ['CLUSTERS','Clusters','clusters']:
            pass

        elif traj in ['NODES','Nodes','nodes']:
            pass
    

    def kinetic_network(self,traj=None,verbose=False):

        if traj in ['CLUSTERS','Clusters','clusters']:
            if type(self.traj_clusters) not in ['numpy.ndarray']:
                self.traj_clusters=array(self.traj_clusters,order="Fortran")
            self.network_clusters=kinetic_network(self.traj_clusters,ranges=[self.traj_clusters.min(),self.traj_clusters.max()],verbose=verbose)
            pass

        elif traj in ['NODES','Nodes','nodes']:
            if type(self.traj_nodes) not in ['numpy.ndarray']:
                self.traj_nodes=array(self.traj_nodes,order="Fortran")
            self.network_nodes=kinetic_network(self.traj_nodes,ranges=[self.traj_nodes.min(),self.traj_nodes.max()],verbose=verbose)
            pass

        else:
            print 'traj= clusters or nodes'
            pass



    def Ganna2012(self,window=None,ksi=0.5,granularity=1.2,bins=20,segment=None,delta_x=None,clusters=True,verbose=False):

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

        if self.num_particles==1:
            self.traj_nodes=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,array([self.traj],order='Fortran'),ksi,window,self.num_particles,self.length)[0]
        else:
            self.traj_nodes=ftrajs.aux.ganna(opt_range,opt,bins,mmn,mmx,delta_x,self.traj,ksi,window,self.num_particles,self.length)


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
