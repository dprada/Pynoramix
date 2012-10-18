.. _getting_started_pyn:

***************
Getting Started
***************

Last developing version
=======================

The whole project is public available at https://github.com/dprada/Pynoramix.git .
There is not a stable version yet, for this reason the use of these libraries it is under your responsability.
New functions or corrections are updated daily. Because of this, keeping an eye on the history of the project is highly recommended:
https://github.com/dprada/Pynoramix/commits/devel



Getting Pynoramix
=================

The last **developing** version can be obtained with any of the following procedures:


*Web*
+++++

The source code can be downloaded in the web page https://github.com/dprada/Pynoramix.git as a zip or tar.gz file.
Both options can be found with the options "zip" or "Downloads":

.. image:: _static/screenshot_github.png


*Git (recommended)*
+++++++++++++++++++

The project can be cloned with git:

.. sourcecode:: bash

   git clone git://github.com/dprada/Pynoramix.git

This former option is recommended. The libraries can be updated easily this way.

.. Todo::
   Adding the project in easy_install or setup.py (http://packages.python.org/an_example_pypi_project/setuptools.html#registering-your-project)
.. Todo::
   Links to raolab.com or GitHub from raolab.



Installing
===========

Pynoramix depends on some packages:

- Python 2.7
- Fortran Compiler (gfortran, intel fortran compiler, ...)
- NumPy
- f2py
- python-dev (to fix the problem with file Python.h)
- liblapack and liblapack-dev (or similar: atlas, blas, mkl, ...)


After solving the dependencies, the Makefile needs to be executed to compile the fortran core of Pynoramix.
This installation script has some variables which can be fullfilled manually:

.. sourcecode:: bash

   F2PY=             # f2py command (f2py,f2py2,...)
   FCOMP=            # fortran compiler command (gfortran, ifort,...)
   FTYPE=            # fortran compiler for f2py (not manually given)
   LAPACK_LIBS=      # lapack libraries (-llapack, -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread, ...)
   FOPTS=            # options of the fortran compiler used (-fast, -checkall, ...)
   FFLAGS=           # additional fortran flags

If these variables are left in blank, they will be detected automatically. 
At this point, and in the directory of Pynoramix, the following command needs to be executed:

.. sourcecode:: bash

   make

If the installation run without troubles, Pynoramix is ready to be used.

.. warning:: Do not forget to add Pynoramix to your python path:
   - export PYTHONPATH=$PYTHONPATH:/path/to/Pynoramix


Being updated
=============

The last modifications can be easily downloaded if you made a git clone.
The command 'git pull' can be executed over the Pynoramix directory to check and obtained the changes.

.. sourcecode:: bash

   git pull

Once this has been done, compiling the changed libraries is mandatory. 
Since the Makefile script detects the changes, running it again is enough:

.. sourcecode:: bash

   make

