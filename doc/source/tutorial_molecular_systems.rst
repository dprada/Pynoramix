
Tutorial Molecular Systems
**************************

First of all, lets load Pynoramix in our script or in a ipython session:

.. sourcecode:: ipython

     In [1]: from pynoramix_beta import *


Some basic notions on python will be assumed along this tutorial. If you just landed here without any idea on python, have a look to the section *First steps on python*.

.. todo:: Make a short tutorial on python, enough to run pynoramix.

----------------------
 

Loading/Writting the topology
=============================

Loading
+++++++

A system can be loaded from a file (pdb,gro) or downloaded from the Protein Data Bank.

.. sourcecode:: ipython

   In [2]: mol_test=molecule(download='2WC2')
   # File saved as 2WC2.pdb
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   # 20  frames/models in traj 0

   In [3]: mol_test=molecule('2WC2.pdb',verbose=False)

   In [4]: mol_test.info()
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   # 20  frames/models in traj 0


.. todo:: Complete the topology of residues and terminals. Apparently
   it works because the option "with_bonds=False" is by default. It needs to be watched.

This way the pdb file has been loaded together with the coordinates
present in the file, 20 different models in the case of 2WC2.  If the
coordinates are going to be loaded from a different file, the topology
can be created adding the option coors=False:

.. sourcecode:: ipython

   In [5]: mol_test2=molecule('2WC2.pdb',coors=False)
   # System created from the file 2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions

   In [6]: mol_test2.info_trajs()
   # No coordinates

   In [7]: mol_test.info_trajs()
   # 20 frames/models in traj 0


Navigating
++++++++++



Writting
++++++++


sdfsdfdsf

----------------------

Loading/Writting a MD trajectory
================================

Once a topology has been created a trajectory can be loaded from
different formats: pdb, gro, xtc, trr, dcd, bin (to be deprecated).

It is recommended the use of dcd files, the file is unformatted and
thereby it is small and easy to handle.

Along this section the different ways to do it will be illustrated.

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb')
   # System created from the file GSGS.pdb :
   # 4723  atoms
   # 1568  residues
   # 3  chains
   # 1560  waters
   # 4  ions
   # 1  frames/models in traj 0

   In [3]: GSGS.delete_traj()
    
   In [4]: GSGS.info_trajs()
   # No coordinates

   In [5]: GSGS.load_traj('GSGS.dcd','ALL')
   # 10 frames/models loaded.

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)
    
   In [3]: GSGS.load_traj('GSGS.dcd',frame='ALL',verbose=False)
    
   In [4]: GSGS.info(); GSGS.info_trajs()
   # System created from the file GSGS.pdb :
   # 4723  atoms
   # 1568  residues
   # 3  chains
   # 1560  waters
   # 4  ions
   # 10  frames/models in traj 0

.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)

   In [3]: GSGS.load_traj('GSGS.dcd')
   # 0 frames/models in traj 0

   In [4]: print GSGS.traj[0].name, GSGS.traj[0].type, GSGS.traj[0].io_opened, GSGS.traj[0].io_end
   GSGS.dcd dcd True False

   In [5]: while 1:
     ....:     GSGS.traj[0].upload_frame()
     ....:     if GSGS.traj[0].io_end: break
     ....: 
   # End of file

   In [6]: GSGS.info_trajs()
   # 10 frames/models in traj 0


.. sourcecode:: ipython

   In [2]: GSGS=molecule('GSGS.pdb',coors=False,verbose=False)

   In [3]: GSGS.load_traj('GSGS.dcd',frame=0)
   # 1 frames/models in traj 0

   In [4]: while GSGS.traj[0].io_opened:
      ...:     print GSGS.traj[0].frame[0].coors[0]
      ...:     GSGS.traj[0].reload_frame()
      ...: 
   [ -7.26851273  -8.12112617  10.57811832]
   [ -5.16595078  -9.8920269   12.24640751]
   [ -6.12880325  -9.20014763  15.28322697]
   [ -4.90646744  -8.31535339  12.97708988]
   [ -5.04781723  -9.68705559  14.15655327]
   [ -5.95707321  -8.45479965  17.51550102]
   [ -4.45994186 -10.63479614  16.19140053]
   [ -6.01659775 -13.60509872  16.98220253]
   [ -4.40946579 -13.10482597  17.12298393]
   [ -5.01924515 -13.77911949  15.64630699]
   # End of file

   In [5]: GSGS.info_trajs()
   # 1 frames/models in traj 0





Converting a trajectory into other format
+++++++++++++++++++++++++++++++++++++++++


How to make an atoms selection
==============================

xxxxxxxxxxxxxx


