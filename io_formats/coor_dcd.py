from frame import *
from numpy import array

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame


##

def open_traj(file_name):
    FFF=libdcd.open(file_name,'r')
    return FFF,0,0


def read_all(file_unit,iopos=None):

    pass
