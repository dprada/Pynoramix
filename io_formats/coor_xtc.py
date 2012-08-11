from frame import *
from numpy import *
from ctypes import *
from numpy.ctypeslib import ndpointer

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, \
exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, \
exdrNOMEM, exdrENDOFFILE, exdrNR = range(13)


try: 
    xdr=cdll.LoadLibrary("libxdrfile.so")
except:
    raise IOError("libxdrfile.so can't be loaded")

def open_traj(file_name):
    FFF=xdr.xdrfile_open(file_name,"r")
    if not FFF: 
        print "Cannot open file: "+file_name
        return FFF,0,1
    else:

        #read natoms
        natoms=c_int()
        r=xdr.read_xtc_natoms(file_name,byref(natoms))
        if r!=exdrOK: raise IOError("Error reading: "+file_name)
        natoms=natoms.value

        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
        xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float), \
           ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]

        return FFF,[natoms],0,0

def read_all(file_unit,iovars=None,iopos=None):

    temp=[]
    while 1:
        temp_frame,pos,end=read_next (file_unit,iovars,iopos)
        if end:
            break
        temp.append(temp_frame)

    return temp,0


def read_next (file_unit,iovars=None,iopos=None):

    step = c_int()
    time = c_float()
    prec = c_float()
    lam = c_float()
    natoms=iovars[0]

    temp_frame=cl_frame()
    temp_frame.coors=empty((natoms,3),dtype=float32)
    temp_frame.box=empty((3,3),float32)

    result = xdr.read_xtc(file_unit,natoms,byref(step),byref(time),temp_frame.box,
                               temp_frame.coors,byref(prec))

    if result!=exdrOK: raise IOError("Error reading xdr file")
    if result==exdrENDOFFILE: 
        print '# End of file'
        return temp_frame,pos,1
     
    temp_frame.precision=prec.value
    temp_frame.time=time.value
    temp_frame.step=step.value
     
    temp_frame.coors=array(temp_frame.coors,order='Fortran')
    temp_frame.coors=10.0*temp_frame.coors
    temp_frame.box=array(temp_frame.box,order='Fortran')
    temp_frame.box=10.0*temp_frame.box

    pos=0
    return temp_frame,pos,0

        
def read_frame(file_unit,iovars=None,iopos=None):

    return None,None,True

def close(file_unit):

    file_unit.close()
    return None,False
    

def write_all():

    pass

def write_frame():

#extern int write_xtc(XDRFILE *xd,
#                       int natoms,int step,float time,
#                       matrix box,rvec *x,float prec);

    pass



