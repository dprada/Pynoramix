Pyno_Back
=========

Real devel version of Pynoramix.

-------------------------------------------

Instructions to compile the fortran subroutines:

With ifort:

f2py -c -m pyn_fort_general pyn_fort_general.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread
f2py -c -m pyn_fort_enm pyn_fort_enm.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread

With gfortran in Goldrake:

f2py2 -c -m pyn_fort_general pyn_fort_general.f90 -llapack
f2py2 -c -m pyn_fort_enm pyn_fort_enm.f90 -llapack

Note: f2py requires python-dev

--------------------------------------------

Provisional Structure:

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
The file pynoramix.py can be exported in an enviroment variable called PYTHONPATH
The part:
import sys
sys.path.append('/home/diego/Projects/Git_Pynoramix/backdoor_Water2/')
... must be included.

