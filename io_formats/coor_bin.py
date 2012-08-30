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

# io_vars[0]: Number of atoms
# io_vars[1]: format frame
# io_vars[2]: box?
# io_vars[3]: pos_header
# io_vars[4]: pos_frame

def open_traj(file_name):

    io_vars=[]
    io_pos=0
    io_err=0

    funit=open(file_name,'rb')
    
    if not funit:
        io_err=1
    else:
        natoms=stc.unpack('i', funit.read(4))[0]
        funit.read(3*4)
        funit.read(3*4)
        funit.read(3*4)
        io_vars.append(natoms)
        io_vars.append(str(3*natoms)+'f')
        io_vars.append([])
        io_pos=funit.tell()
        io_vars.append(io_pos)
        io_vars.append((13+natoms*3)*4)

    return funit,io_vars,io_pos,io_err

def read_all(file_unit,io_vars=None,io_pos=None):

    temp=[]
    while 1:
        temp_frame,io_pos,io_err,io_end=read_next(file_unit,io_vars,io_pos)
        if io_end or io_err:
            break
        temp.append(temp_frame)

    return temp,io_err,io_end   # io_file,io_err,io_end

def read_next(file_unit,io_vars=None,io_pos=None):

    temp_frame,io_pos,io_err,io_end=read_aux(file_unit,io_vars,io_pos)
    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def read_frame(file_unit,frame,io_vars=None,io_pos=None):

    io_pos=io_vars[3]+frame*io_vars[4]
    temp_frame,io_pos,io_err,io_end=read_aux(file_unit,io_vars,io_pos)
    return temp_frame,io_pos,io_err,io_end

def read_aux(file_unit,io_vars=None,io_pos=None):

    io_err=0
    io_end=0
    temp_frame=cl_frame()

    file_unit.seek(io_pos)

    test_byte=file_unit.read(1)
    if test_byte=='': 
        io_end=1
    else:
        file_unit.seek(io_pos)

        natoms=stc.unpack('i', file_unit.read(4))[0]
        natoms3=natoms*3
        temp_frame.step=stc.unpack('i', file_unit.read(4))[0]
        temp_frame.time=stc.unpack('f', file_unit.read(4))[0]
        temp_frame.box[0][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]
        temp_frame.box[1][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]
        temp_frame.box[2][:]=stc.unpack('3f', file_unit.read(3*4))[0:3]                                                 
        temp=stc.unpack(io_vars[1], file_unit.read(natoms3*4))[0:natoms3]
        for ii in range(0,natoms3,3):
            aux=temp[ii:ii+3]
            temp_frame.coors.append(aux)
        temp_frame.box=10.0*array(temp_frame.box,order='Fortran')
        temp_frame.coors=10.0*array(temp_frame.coors,order='Fortran')
        temp_frame.precision=stc.unpack('f',file_unit.read(4))[0]             # Precision of the trajectory

        io_pos=file_unit.tell()

    return temp_frame,io_pos,io_err,io_end  # frame,io_pos,io_err,io_end

def close_traj(file_unit):

    file_unit.close()
    io_err=0
    return io_err

def write_all():
    print '# Not Supported'
    pass

def write_frame():
    print '# Not Supported'
    pass

