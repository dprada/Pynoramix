Tutorial Kinetic Analysis
*************************

First of all, lets load Pynoramix in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from pynoramix_beta import *


Some basic notions on python will be assumed along this tutorial. If you just landed here without any idea on python, have a look to the section *First steps on python*.

.. todo:: Make a short tutorial on python, enough to run pynoramix.

----------------------

1D Trajectories
===============

Given a trajectory, we will temporary work with 1D traj, these are
some functions for its analysis.

Anna's method
+++++++++++++

It should be explained.

.. sourcecode:: ipython

   In [2]: traj_ganna=ganna(traj1D,window=tw,ksi=0.5,delta_x=0.20,segment=[-21.0,21.0])
    
   In [3]: net=kinetic_network(traj_ganna,ranges=[min(traj_ganna),max(traj_ganna)],verbose=False)
    
   In [4]: net_s=net.symmetrize(verbose=False)
    
   In [5]: net_s.mcl(granularity=1.2,pruning=True,verbose=False)
    
   In [6]: traj_mcl_clusts=[[] for ii in range(net_s.num_clusters)]
    
   In [7]: for ii in range(len(traj_ganna)):
     ....:     nn=traj_ganna[ii]
     ....:     cluster=net_s.node[nn].cluster
     ....:     traj_mcl_clusts[cluster].append(traj1D[tw+ii])
     ....: 


Rao's method
++++++++++++

.. sourcecode:: ipython

   In [2]: traj_rao=rao(traj1D,window=25,separators=[-6.00,-2.40])
    
   In [3]: net,traj_nodes=kinetic_network(traj_rao,ranges=[[0,51],[0,51],[0,51]],traj_out=True,verbose=False)
    
   In [4]: len(traj_nodes[0])
   Out[4]: 9998951
    
   In [5]: net.node[0].label
   Out[5]: '[ 0  0 51]'
    
   In [6]: net.info()
   # Network:
   # 754 nodes
   # 2618 links out
   # 9998950.0 total weight nodes

   In [7]: net_s=net.symmetrize(verbose=False)

   In [8]: net_s.mcl(granularity=1.2,pruning=True,verbose=True)

   In [9]: traj_mcl_clusts=[[] for ii in range(net_s.num_clusters)]
    
   In [10]: for ii in range(len(traj_nodes[0])):
      ....:       cluster=net_s.node[traj_nodes[0][ii]].cluster
      ....:       traj_mcl_clusts[cluster].append(traj1D[tw+ii])
      ....: 





