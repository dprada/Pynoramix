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

Given a 1D trajectory , its kinetic analysis can be done taking advantage of the class *kinetic_1D_analysis*. 
The following subsections are a description of how to deal with this class and the functions included.

.. seealso:: The section XXX for the details of this function.


Loading a trajectory
++++++++++++++++++++

Given a file, :download:`traj.oup
<../tutorials/kinetic_1D_analysis/traj.oup.gz>`, with the 1D
trajectory in its second column, the class *kinetic_1D_analysis* can
be initialized:

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

First order kinetic model
+++++++++++++++++++++++++

At first glance the former distribution can be decompose in three
peaks, called macro-states since we know that this histogram is coming
from a dynamical system. Actually, the trajectory corresponds to the
projection onto an arbitrary coordinate of a langevin dynamics of a
particle embeded in a 2D potential.

.. figure:: ../tutorials/kinetic_1D_analysis/histo_color.png
   :align: center
   :scale: 70 %

According to this first superficial analysis, the particle moves
through a landscape with 3 basins (green, orange, blue) separated by 2
barriers (green-orange and orange-blue). Thereby, the particle goes
from green to blue region bymeans of a necessary visit to the region orange. 

A first order kinetic model can be computed with these three
states. The model should be written as a master equation where the
transition probabilities, or rates, between macrostates interplay to
mimic the macroscopical dynamical behavior observed in the experiment.

These transition probabilities are computed from the following trajectory:

.. sourcecode:: ipython

   In [6]: for xx in kin_test.traj[:]:
      ...:     if xx <=-5.5:
      ...:          kin_test.traj_clusters.append(0)   # Green region
      ...:     elif xx <=-2.3:
      ...:          kin_test.traj_clusters.append(1)   # Orange region
      ...:     else:
      ...:          kin_test.traj_clusters.append(2)   # Blue region
      ...: 

Notice that the term 'cluster' was used: a set of time steps belonging
to the same macrostate. This way the cluster labeled '0' corresponds
to region green, '1' to the orange one and '2' to the blue.  The
transition probabilities can be computed building the corresponding
kinetic network, or conformational space network. The trajectory is
mapped into a network where weighted nodes and directed links accounts to
the occupation of the clusters and the number of transitions among them.

.. sourcecode:: ipython

   In [7]: kin_test.kinetic_network(traj='clusters',verbose=True)
   # Network:
   # 3 nodes
   # 9 links out
   # 999900.0 total weight nodes

Before getting the transision probabilities, or rates, detailed
balance condition is impossed symmetrising the network.

.. sourcecode:: ipython

   In [8]: kin_test.network_clusters.symmetrize(new=False,verbose=False)

   In [9]: kin_test.network_clusters.info()
   # Network:
   # 3 nodes
   # 9 links out
   # 1999800.0 total weight nodes

.. seealso:: The section XXX for the details of this function.


The network we have now is nothing but a kinetic model respresented as a graph. This way the master equations:

.. math::

   \frac{dP_{0}}{dt} = P_{01}P_{1} + P_{02}P_{2} - P_{10}P_{0} - P_{20}P_{0},

\

.. math::

   \frac{dP_{1}}{dt} = P_{10}P_{0} + P_{12}P_{2} - P_{01}P_{1} - P_{21}P_{1},

\

.. math::

   \frac{dP_{2}}{dt} = P_{20}P_{0} + P_{21}P_{1} - P_{02}P_{2} - P_{12}P_{2},


can be written with the transition probabilities:

.. sourcecode:: ipython

   In [10]: for ii in range(3):
      ....:    ww=kin_test.network_clusters.node[ii].weight
      ....:        for jj in range(3):
      ....:            if ii!=jj :
      ....:                print 'P'+str(jj)+str(ii)+'=', kin_test.network_clusters.node[ii].link[jj]/ww
      ....: 
   P10= 0.00134372002214
   P20= 0.000153059499211
   P01= 0.00190455529265
   P21= 0.0749385749386
   P02= 8.25724981733e-05
   P12= 0.0285230372469

Where detailed balance condition, :math:`P_{ji}P^{s}_{i}=P_{ij}P^{s}_{j}`, is
fullfilled by construction. Thereby, the stationary solution of this
master equation is given by the occupation probabilities given found
as weights of nodes.

.. sourcecode:: ipython

   In [11]: ww=kin_test.network_clusters.weight

   In [12]: for ii in range(3):
      ....:     print 'P'+str(ii)+'^s =', kin_test.network_clusters.node[ii].weight/ww
      ....: 
   P0^s = 0.28096459646
   P1^s = 0.198228822882
   P2^s = 0.520806580658



We can now compare our kinetic model with the behavior observed in the
original trajectory.  This can be done attending to magnitudes as the
*mean dwell/life time* and its distribution or the *mean first passage
time* and its distribution.

Mean life time or life time distribution
++++++++++++++++++++++++++++++++++++++++

Up to here, in this tutorial we have built a kinetic model from a 1D
trajectory. Now the accuracy of the model can be checked with
observables as the life time of its macro-states.

First of all, the life time distribution of the 3 macro-states can be
computed from our clusters trajectory.

.. sourcecode:: ipython

   In [13]: ltx_1,lty_1=kin_test.life_time(traj='clusters',state=1,norm=True,verbose=True)
   # Mean life time: 16.3446993094 frames.

   In [14]: ltx_2,lty_2=kin_test.life_time(traj='clusters',state=2,norm=True,verbose=True)
   # Mean life time: 34.0213561056 frames.

.. Note:: This analysis was applied not to the trajectory found in the
   file 'traj.oup' but to a trajectory 10 times larges.

.. Seealso:: Function for further details.

The distributions can be fitted to exponential functions:
:math:`\exp{(At+B)}`. The next plot shows the fitting together with
the distributions, where approximately A=-0.008 and B=-9.763 for
cluster 1 (orange), and A=-0.026 and B=-7.311 for cluster 2 (blue).

.. figure:: ../tutorials/kinetic_1D_analysis/lf_dist_bad.png
   :align: center
   :scale: 70 %

We observe how the long time behavior fits to an exponential
distribution but **many non expected short life events appear**.

Now, we can compare with the life time distributions obtained by the
kinetic model although from the previous plot we can guess that our
clusters are not well defined. To this aim a random walker over the
kinetic network can be run according to the transition probabilities.

.. sourcecode:: ipython

   In [15]: bw_traj=kin_test.network_clusters.brownian_walker(origin=0,length=10000000)

   In [16]: bw=kinetic_1D_analysis(bw_traj)
   # Trajectory loaded: 10000001 time steps

   In [17]: ltbw_x_1,ltbw_y_1=bw.life_time(state=1,norm=True,verbose=True)
   # Mean life time: 16.2835161764 frames.

   In [18]: ltbw_x_2,ltbw_y_2=bw.life_time(state=2,norm=True,verbose=True)
   # Mean life time: 34.1362159764 frames.

This time the distributions and their exponential fits, A=-0.068 and B=-5.504 for
cluster 1 (orange) and A=-0.029 and B=-7.034 for cluster 2 (blue), shows a different behavior:

.. figure:: ../tutorials/kinetic_1D_analysis/lf_dist_model_bad.png
   :align: center
   :scale: 70 %


We can see how, **although the mean life time (in [13],[14] and
[17],[18]) are similar, the relaxations times are different from
those previously computed**. Cluster 1 (orange) has a faster
relaxation than Cluster 2 (blue). It is up to the reader trying to
understand the origin of these different behaviors, even with the same
*mlt* (check section XXX to find some help).


These observations point to the possibility of having a bad kinetic
description of our system, but before showing how to do it better,
lets check another kinetic magnitude: the *first passage time distribution*.



Mean first passage time or fpt distribution
+++++++++++++++++++++++++++++++++++++++++++

Following the same strategy as in the previous section, we can compare
some distributions computed over the original clusters trajectory and
over its kinetic model. The first passage time distribution can be
computed to a certain state from a given state or from any one.

The original clusters trajectory shows the following distributions for
this three different cases: from any cluster to cluster 2 (blue), from
cluster 0 (green) to cluster 2 and from cluster 1 (orange) to cluster 2.

.. sourcecode:: ipython

   In [19]: fptx2,fpty2=kin_test.first_passage_time(traj='clusters',to_state=2,norm=False,verbose=True)
   # Mean first passage time: 700.403755073 frames.

   In [20]: fptx02,fpty02=kin_test.first_passage_time(traj='clusters',from_state=0,to_state=2,norm=False,verbose=True)
   # Mean first passage time: 1175.03303483 frames.

   In [21]: fptx12,fpty12=kin_test.first_passage_time(traj='clusters',from_state=1,to_state=2,norm=False,verbose=True)
   # Mean first passage time: 124.958081567 frames.


.. figure:: ../tutorials/kinetic_1D_analysis/fpt_bad_to2.png
   :align: center
   :scale: 70 %

There are two things not expected in this plot. The **short time
behavior looks anomalous** and both distributions '1 to 2' and '0 to
2' should be bi-exponential. Since we made the first guest of having the
particle going from **0 to 2** with a necessary stay in cluster 1, this
distribution **should be bi-exponential**.


Now, we can compare with the first passage time distributions obtained
by the kinetic model. As it was done with the life time distributions,
a random walker over the kinetic network can be run according to the
transition probabilities.

.. sourcecode:: ipython

   In [22]: bw_traj=kin_test.network_clusters.brownian_walker(origin=0,length=9999000)

   In [23]: bw=kinetic_1D_analysis(bw_traj,verbose=False)

   In [24]: fptbwx2, fptbwy2= bw.first_passage_time(to_state=2,norm=False,verbose=True)
   # Mean first passage time: 430.700895225 frames.

   In [25]: fptbwx02,fptbwy02=bw.first_passage_time(from_state=0,to_state=2,norm=False,verbose=True)	
   # Mean first passage time: 751.240082388 frames.

   In [26]: fptbwx12,fptbwy12=bw.first_passage_time(from_state=1,to_state=2,norm=False,verbose=True)
   # Mean first passage time: 36.4134420305 frames.


.. figure:: ../tutorials/kinetic_1D_analysis/fpt_bad_to2_model.png
   :align: center
   :scale: 70 %


This time we atleast observe a bi-exponential behavior in the '1 to 2'
distribution in agreement with our first guess. But if we compare the
relaxation times, **by simple inspection we conclude that distributions
are different** from those obtained before from the clusters trajectory,
**even the mean first passage times are different (in [19-21] and
[24-26])**.

These analysis reinforce the idea of having a wrong analysis, the way
we have define the macro-states or cluster is not kinetically
consistent.  But before improving this analysis, let see one kinetic
observable more: what we have define as *first committed passage time*
(from the magnitude: committor probabilities)

It is up to the reader trying to understand the origin of these
different behaviors (check section XXX to find some help).



Mean first committed passage time or fcpt distribution
++++++++++++++++++++++++++++++++++++++++++++++++++++++

In Physics of stochastic processes the *committor probability*,
*C_{a,b}*, is the probability of going from *b* to *a* without passing
through any other state. Based on this idea, we can compute the first
passage time from a state *b* to *a* with a requirement: passing or
not by other states.

In this tutorial we defined 3 macro-states or clusters from a 1D
trajectory (see above). And in our naive first analysis we made an
observation: to go from cluster 0 (green) to cluster 2 (blue), a visit
to cluster 1 (orange) is needed. For the sake of completeness lets
study the first passage time of going from 0 to 2 with and without
passing by cluster 1.

fcpt_without1=kin_test.first_committed_passage_time(traj='clusters',states=[0,2],verbose=True)
fcpt_with1=   kin_test.first_committed_passage_time(traj='clusters',states=[0,1,2],commitment=[True,False,True],verbose=True)
fcpt_without1=kin_test.first_committed_passage_time(traj='clusters',states=[0,1,2],commitment=[True,True,True],verbose=True)



.. figure:: ../tutorials/kinetic_1D_analysis/traj123.png
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





