Ion Solvation
=============

This is the analysis of a water box with one ion.

Loading the system and trajectory.
**********************************

The whole system is loaded for its analysis (I should rewritte the
analysis reading frame by frame).

.. sourcecode:: ipython

   In [1]: from pynoramix_beta2 import *
    
   In [2]: import pylab as pylab
    
   In [3]: ion=molecule('run_ion.gro',coors=False,verbose=False) 
    
   In [4]: ion.load_traj('traj.dcd',frame='ALL',verbose=False)
    
   In [5]: list1=ion.selection('atom.resid.type Ion')
    
   In [6]: list2=ion.selection('atom.name OW')

   In [7]: list3=ion.selection('atom.type H')

   In [8]: ion.info(); ion.info_trajs()
   # System created from the file run_ion.gro :
   # 4093  atoms
   # 1024  residues
   # 0  chains
   # 1023  waters
   # 1  ions
   # 52501 frames/models in traj 0

   In [9]: print '#',len(list1), 'atoms selected in list1'
   # 1 atoms selected in list1
    
   In [10]: print '#',len(list2), 'atoms selected in list2'
   # 1023 atoms selected in list2

   In [11]: print '#',len(list3), 'atoms selected in list3'
   # 2046 atoms selected in list3


Radial Distribution Function
****************************

RDFs: Ion-Oxygens and Ion-Hydrogens.

.. sourcecode:: ipython

   In [8]: rdf_x_o,rdf_y_o=ion.rdf(setA=list1,setB=list2,bins=1000,segment=[0.0,30.0])
    
   In [9]: rdf_x_h,rdf_y_h=ion.rdf(setA=list1,setB=list3,bins=1000,segment=[0.0,30.0])

   In [10]: pylab.plot(rdf_x_o,rdf_y_o,'r-',rdf_x_h,rdf_y_h,'b-')

The cut-offs from 1st shell and 2nd shell for O are: 3.6 and 5.75
The cut-offs from 1st shell and 2nd shell for H are: 4.25 and 6.50


Definition of bond Ion-Water (oxygen and or hydrogen)
*****************************************************

.. sourcecode:: ipython

   In [10]: distances=ion.distance(list2,list1)

   In [11]: total_net=network(directed=True,kinetic=True,verbose=False)
    
   In [12]: for ii in range(1023):
      ....:     traj_states=rao(distances[:,ii,0],window=25,separators=[3.600,5.750])
      ....:     net,traj_nodes=kinetic_network(traj_states,ranges=[[0,51],[0,51],[0,51]],traj_out=True,verbose=False)
      ....:     total_traj_nodes.append(traj_nodes)
      ....:     total_net.merge_net(net,verbose=False)
      ....: 

   In [13]: total_net_s=total_net.symmetrize(verbose=False)     # I should change this
    
   In [14]: total_net_s.mcl(granularity=1.3,pruning=True,verbose=True)
   # Number of clusters:  3

   In [15]: traj_mcl_clusts=[[] for ii in range(total_net_s.num_clusters)]
    
   In [16]: tw=25
    
   In [17]: for jj in range(1023):
      ....:         for ii in range(len(total_traj_nodes[jj])):
      ....:           cluster=total_net_s.node[total_traj_nodes[jj][ii]].cluster
      ....:       traj_mcl_clusts[cluster].append(distances[tw+ii,jj,0])
      ....: 

