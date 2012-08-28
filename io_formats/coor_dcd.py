from frame import *
from numpy import array
import libdcdfile as libdcd

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame


##

# FFF: unit number for Fortran

def open_traj(file_name):
    funit=libdcd.open(len(file_name),str(file_name),'r')

    if not funit:
        print "Cannot open file: "+file_name
        return funit,0,1
    else:
        return funit,0,0


def read_all(file_unit,iopos=None):

    temp=[]
    while 1:
        temp_frame,pos,end=read_next(file_unit,iovars,iopos)
        if end:
            break
        temp.append(temp_frame)

    return temp,0

def read_next(file_unit,iovars=None,iopos=None):

    return None,None,1

def read_frame(file_unit,iopos=None):

    return None,None,1

def close_traj(file_unit):

    file_unit.close()
    return None,0
    

def write_all():

    pass

def write_frame():

    pass
