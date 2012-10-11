Tutorial Kinetic Analysis
*************************

First of all, lets load Pynoramix in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from pynoramix_beta import *


Some basic notions on python will be assumed along this tutorial. 


.. If you just landed here without any idea on python, have a look to the section *First steps on python*.
.. coment todo:: Make a short tutorial on python, enough to run pynoramix.

----------------------

1D Trajectories
===============

Given a trajectory 1D trajectory, its kinetic analysis can be done taking advantage of the class *kinetic_1D_analysis*. 
The following subsections are a description of how to deal with this class and the functions included.

.. seealso:: The section XXX for the details of this function.


Loading a trajectory
++++++++++++++++++++

Given a file named "traj.oup" with the 1D trajectory in its second column, the class *kinetic_1D_analysis* can be initialized:

.. sourcecode:: ipython

   In [2]: kin_test=kinetic_1D_analysis('traj.oup',column=1)
   # Trajectory loaded: 999901 time steps

Now, the trajectory is stored in the object kin_test:

.. sourcecode:: ipython

   In [3]: print 'Column', kin_test.file_column, 'in file', kin_test.file_name
   Column 1 in file traj.oup

   In [4]: print kin_test.num_particles, 'particles  with', kin_test.dimensions, 'dimension', ':', kin_test.length, 'time steps.'
   1 particles with 1 dimension: 999901 time steps.

   In [5]: print kin_test.traj[:]
   [-0.91936072  0.74886578 -1.07525923 ..., -9.09667346 -8.40884264 -8.83918787]

Since the trajectory is a numpy.ndarray, we can take adventage of its numpy intrinsic attributes and functions:

.. sourcecode:: ipython

   In [6]: print 'Min.:',kin_test.traj.min(),'   Mean:', kin_test.traj.mean() ,'   Max.:', kin_test.traj.max()
   Min.: -11.4344662381    Mean: -3.40808122446    Max.: 2.78106514618

The histogram can be obtained with the following command:

.. sourcecode:: ipython

   In [7]: hxx,hyy=kin_test.histogram(delta_x=0.20,segment=[-12.0,4.0],norm=False)

   In [8]: import pylab as pylab

   In [9]: pylab.plot(hxx,hyy,'b-')
   Out[10]: [<matplotlib.lines.Line2D object at 0x43277d0>]

   In [11]: pylab.show()

.. figure:: ../tutorials/kinetic_1D_analysis/histo_1D.png
   :align: center
   :scale: 70 %

Accurate kinetic decomposition
++++++++++++++++++++++++++++++


Ganna2012
.........

The following method implements the analysis propossed by
G. Berezovska et al. [CITE] to unveil the conformational macrostates
kinetically "well defined" and the underlying accurate first order
kinetic model.

.. sourcecode:: ipython

   In [3]: kin_test.Ganna2012(window=25,granularity=1.2,bins=20,verbose=True)
   # Network:
   # 1260 nodes
   # 49972 links out
   # 1999700.0 total weight nodes
   # Number of clusters:  3

.. seealso:: The section XXX for the details of this function.




.. Warning::

   Please cite the following reference if the method is used for a scientific publication: XXXXXXX

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





