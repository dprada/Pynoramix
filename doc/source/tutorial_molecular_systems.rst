
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

   

.. todo:: Complete the topology of residues and terminals. Apparently
   it works without the option "with_bonds=False" but it needs to be watched.

Navigating
++++++++++

Writting
++++++++


sdfsdfdsf

----------------------

Loading/Writting a MD trajectory
================================

xxxxxxxxxxxxxxxxxxxx

Converting a trajectory into other format
+++++++++++++++++++++++++++++++++++++++++


How to make an atoms selection
==============================

xxxxxxxxxxxxxx


