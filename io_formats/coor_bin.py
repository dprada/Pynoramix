from frame import *
from numpy import array
import struct as stc #To read binary files

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

### This format should be deprected.

#if self.type in ['bin','gro','xtc']:
#    self.coors=10.0*self.coors
#    self.box=10.0*self.box


def open_traj(file_name):
    FFF=open(file_name,'rb')
    # Header
    natoms=stc.unpack('i', FFF.read(4))[0]
    FFF.read(3*4)
    FFF.read(3*4)
    FFF.read(3*4)
    format=str(3*natoms)+'f'
    return FFF,[natoms,format],0,0

def read_all(file_unit,iovars=None,iopos=None):

    temp=[]
    pos=0
    file_unit.seek(pos)
    natoms=iovars[0]
    format=iovars[1]
    
    # Header:
    natoms=stc.unpack('i', file_unit.read(4))[0]
    file_unit.read(3*4)
    file_unit.read(3*4)
    file_unit.read(3*4)
    #self.box[0][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]   
    #self.box[1][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]
    #self.box[2][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]         
             
    #if frame!= None: 
    #    self.frame=frame            # Going to the choosen frame
    #    bytes_frame=(13+N_A*3)*4       # Bytes per frame
    #    self.file_traj.seek(bytes_frame*self.frame,1)     # Jumping to the choosen frame
    #else:
    #    self.frame+=1

    pos=file_unit.tell()
    end=1
    while END:
        file_unit.seek(pos)
        temp_frame=cl_frame()
        natoms=stc.unpack('i', file_unit.read(4))[0]                              # N_A  = number of atoms
        temp_frame.step=stc.unpack('i', file_unit.read(4))[0]
        temp_frame.time=stc.unpack('f', file_unit.read(4))[0]
        temp_frame.box[0][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        temp_frame.box[1][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]    # for the size of the box vg.box        
        temp_frame.box[2][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]                                                 
        temp=stc.unpack(format, file_unit.read(natoms*4*3))[0:natoms*3]
        for ii in range(0,3*natoms,3):
            aux=temp[ii:ii+3]
        temp_frame.coors=aux
        temp_frame.box[0][0]=10.0*float(line[0])
        temp_frame.box[1][1]=10.0*float(line[1])
        temp_frame.box[2][2]=10.0*float(line[2])
        temp_frame.coors=array(temp_frame.coors,order='Fortran')
        temp_frame.coors=10.0*temp_frame.coors
        temp_frame.precision=stc.unpack('f',file_unit.read(4))[0]             # Precision of the trajectory
        temp.append(temp_frame)
        pos=file_unit.tell()
        test_byte=file_unit.read(1)
        if byte=='': end=0

    return temp,0

def read_next(file_unit,iovars=None,iopos=None):

    natoms=iovars[0]
    format=iovars[1]
    temp_frame=cl_frame()
    natoms=stc.unpack('i', file_unit.read(4))[0]                              # N_A  = number of atoms
    temp_frame.step=stc.unpack('i', file_unit.read(4))[0]
    temp_frame.time=stc.unpack('f', file_unit.read(4))[0]
    temp_frame.box[0][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
    temp_frame.box[1][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]    # for the size of the box vg.box        
    temp_frame.box[2][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]                                                 
    temp=stc.unpack(format, file_unit.read(natoms*4*3))[0:natoms*3]
    for ii in range(0,3*natoms,3):
        aux=temp[ii:ii+3]
        temp_frame.coors.append(aux)
    temp_frame.box=10.0*temp_frame.box
    temp_frame.coors=10.0*array(temp_frame.coors,order='Fortran')
    temp_frame.precision=stc.unpack('f',file_unit.read(4))[0]             # Precision of the trajectory
    pos=file_unit.tell()

    return temp_frame,pos,0

def read_frame(file_unit,iopos=None):

    return None,None,True

def close(file_unit):

    file_unit.close()
    return None,False
    

def write_all():

    pass

def write_frame():

    pass

