Pynoramix
=========

There is not a stable version yet, for this reason the use of these
libraries it is under your responsability.  New functions or
corrections are updated daily. Because of this, keeping an eye on the
history of the project is highly recommended.

Further information in: http://dprada.github.com/Pynoramix/

-------------------------------------------

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


- F2PY=             # f2py command (f2py,f2py2,...)
- FCOMP=            # fortran compiler command (gfortran, ifort,...)
- FTYPE=            # fortran compiler for f2py (not manually given)
- LAPACK_LIBS=      # lapack libraries (-llapack, -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread, ...)
- FOPTS=            # options of the fortran compiler used (-fast, -checkall, ...)
- FFLAGS=           # additional fortran flags


If these variables are left in blank, they will be detected automatically. 
At this point, and in the directory of Pynoramix, the following command needs to be executed:

$ make

If the installation run without troubles, Pynoramix is ready to be used.

Do not forget to add Pynoramix to your python path!!

$ export PYTHONPATH=$PYTHONPATH:/path/to/Pynoramix


--------------------------------------------

Being updated
=============

The last modifications can be easily downloaded if you made a git clone.
The command 'git pull' can be executed over the Pynoramix directory to check and obtained the changes.

$ git pull

Once this has been done, compiling the changed libraries is mandatory. 
Since the Makefile script detects the changes, running it again is enough:

$ make

--------------------------------------------

Structure of code:

>> pynoramix.py   (main file with heads, variables and includes)

   >> pyn_cl_set        (classes -topology- and functions of general purposse)
   >> pyn_cl_coors      (class coordinates)
   >> pyn_fort_general  (fortran subroutines of general purposse)

   >> pyn_cl_enm        (module with elastic network models)
   >> pyn_fort_enm      (specific fortran subroutines for enm analysis)

   >> pyn_cl_water      (module with specific water analysis) -To be added-
   >> pyn_fort_water    (specific fortran subroutines for water analysis) -To be added-

   >> pyn_cl_net        (module with complex network analysis) 
   >> pyn_fort_net      (specific fortran subroutines for network analysis) 

---------------------------------------------