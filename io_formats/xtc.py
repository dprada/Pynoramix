#if self.type in ['bin','gro','xtc']:
#    self.coors=10.0*self.coors
#    self.box=10.0*self.box


###################################################
######## XDRFILES (XTC and TRR)
###################################################
        
from ctypes import *
import os.path
from numpy.ctypeslib import ndpointer

class frame:
    #variables
    #x: rvec*natoms / numpy array if installed
    #box DIM*DIM
    #step 
    #time 
    #prec 
    #lam: lambda

    def __init__(self,n):
        #create vector for x
        self.x=empty((n,3),dtype=float32)
        self.box = empty((3,3),float32)


class xdrfile(labels_file_traj):
    exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC, exdrNOMEM, exdrENDOFFILE, exdrNR = range(13)

    #
    def __init__(self,fn,ft="Auto"):

        self.restart()
        self.name=fn
        self.status='OPENED'

        if ft=="Auto":
          ft = os.path.splitext(fn)[1][1:]
        self.mode=ft

        #load libxdrfil
        try: 
          self.xdr=cdll.LoadLibrary("libxdrfile.so")
        except:
          raise IOError("libxdrfile.so can't be loaded")
          
        #open file
        self.xd = self.xdr.xdrfile_open(fn,"r")
        if not self.xd: raise IOError("Cannot open file: '%s'"%fn)
        
        #read natoms
        natoms=c_int()
        if self.mode == 'trr':
            r=self.xdr.read_trr_natoms(fn,byref(natoms))
        else:
            r=self.xdr.read_xtc_natoms(fn,byref(natoms))
        if r!=self.exdrOK: raise IOError("Error reading: '%s'"%fn)
        self.natoms=natoms.value
        
        #for NumPy define argtypes - ndpointer is not automatically converted to POINTER(c_float)
        #alternative of ctypes.data_as(POINTER(c_float)) requires two version for numpy and c_float array
        self.xdr.read_xtc.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),
           ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float)]
        self.xdr.read_trr.argtypes=[c_int,c_int,POINTER(c_int),POINTER(c_float),POINTER(c_float),
           ndpointer(ndim=2,dtype=float32),ndpointer(ndim=2,dtype=float32),POINTER(c_float),POINTER(c_float)]

    def __iter__(self):
        f = frame(self.natoms)
        #temporary c_type variables (frame variables are python type)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        while 1:
            #read next frame
            if self.mode=='xtc':
                result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,
                        f.x,byref(prec))
                f.prec=prec.value
            else:
                result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),
                        f.box,f.x,None,None) #TODO: make v,f possible
                f.lam=lam.value
                
            #check return value
            if result==self.exdrENDOFFILE: break
            if result==self.exdrINT and self.mode=='trr': 
              break  #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
            if result!=self.exdrOK: raise IOError("Error reading xdr file")
            
            #convert c_type to python 
            f.step=step.value
            f.time=time.value
            yield f

    def upload_frame(self):
        f = frame(self.natoms)
        step = c_int()
        time = c_float()
        prec = c_float()
        lam = c_float()
        if self.mode=='xtc':
            result = self.xdr.read_xtc(self.xd,self.natoms,byref(step),byref(time),f.box,
                     f.x,byref(prec))
            f.prec=prec.value
        else:
            result = self.xdr.read_trr(self.xd,self.natoms,byref(step),byref(time),byref(lam),
                     f.box,f.x,None,None) #TODO: make v,f possible
            f.lam=lam.value
        #check return value
        if result==self.exdrENDOFFILE: 
            print '# End of file',self.name
            self.status='END'
            return
        if result==self.exdrINT and self.mode=='trr': 
            print 'mmm...'
            return 0
            #TODO: dirty hack. read_trr return exdrINT not exdrENDOFFILE
        if result!=self.exdrOK: raise IOError("Error reading xdr file")
            
        #convert c_type to python 
        f.step=step.value
        f.time=time.value
        return f
