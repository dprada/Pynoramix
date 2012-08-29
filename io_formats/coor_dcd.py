from frame import *
from numpy import array
import libdcdfile as libdcd

# io_vars[0]: Number of frames in the file
# io_vars[1]: Number of previous integration steps
# io_vars[2]: Frequency (integration steps) to save this file
# io_vars[3]: Number of integration steps in the run to create this file
# io_vars[4]: Frequency of coordinate saving
# io_vars[5]:
# io_vars[6]:
# io_vars[7]: Number of degrees of freedom during the run
# io_vars[8]: Number of fixed atoms
# io_vars[9]: Timestep in AKMA-units. Bit-copy from the 32-bit real number
# io_vars[10]: 1 if crystal lattice information is present in the frames
# io_vars[11]: 1 if this is a 4D trajectory
# io_vars[12]: 1 if fluctuating charges are present
# io_vars[13]: 1 if trajectory is the result of merge without consistency checks
# io_vars[14]:
# io_vars[15]:
# io_vars[16]:
# io_vars[17]:
# io_vars[18]:
# io_vars[19]: CHARMM version number
#---
# io_vars[20]: Number of atoms
# io_vars[21]: delta_t


# open
# > io_file,io_vars,io_pos,io_err

# close
# read_all
# > [frames],io_err

# read_next
# > temp_frame,io_pos,io_err

# read_frame
# write_all
# write_frame


##

# FFF: unit number for Fortran

def open_traj(file_name):

    io_vars=[]
    io_pos=0
    io_err=False

    funit,o_vars,o_natom,o_delta_t,io_pos=libdcd.open_read(len(file_name),str(file_name))

    if not funit:
        io_err=True
    else:
        for ii in o_vars:
            io_vars.append(ii)
        io_vars.append(o_natom)
        io_vars.append(o_delta_t)

    return funit,io_vars,io_pos,io_err  # io_file,io_vars,io_pos,io_err


def read_all(file_unit,io_vars=None,io_pos=None):

    temp=[]
    while 1:
        temp_frame,io_pos,io_err,io_end=read_next(file_unit,io_vars,io_pos)
        if io_end or io_err:
            break
        temp.append(temp_frame)

    return temp,io_err,io_end   # io_file,io_err,io_end

def read_next(file_unit,io_vars=None,io_pos=None):

    temp_frame=cl_frame()
    #temp_frame.coors=empty((natoms,3),dtype=float32)
    #temp_frame.box=empty((3,3),float32)

    io_pos,temp_frame.box,temp_frame.coors,io_err,io_end=libdcd.read(file_unit,io_vars[20],io_vars[10],io_pos)

    return temp_frame,io_pos,io_err,io_end  # io_file,io_pos,io_err,io_end

def read_frame(file_unit,io_pos=None):

    return None,None,1

def close_traj(file_unit):

    io_err=libdcd.close(file_unit)
    return io_err   #io_err=0 good
    

def write_all():

    pass

def write_frame():

    pass
