
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

   In [2]: mol_test=molecule(download='2WC2',coors=False,with_bonds=False)
   # System created from the file  2WC2.pdb :
   # 6704  atoms
   # 418  residues
   # 2  chains
   # 0  waters
   # 0  ions
   

.. todo:: Complete the topology of residues and terminals and remove the option "with_bonds=False".

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


