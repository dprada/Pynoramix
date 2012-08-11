from frame import *
from numpy import array

# open
# close
# read_all
# read_next
# read_frame
# write_all
# write_frame

def open_traj(file_name):
    FFF=open(file_name,'r')
    return FFF,[],0,0

def read_all(file_unit,iovars=None,iopos=None):

    temp=[]
    model=0
    pos=0
    end=1
    while end:
        temp_frame=cl_frame()
        file_unit.seek(pos)
        line=file_unit.readline()
        line=file_unit.readline()
        num_atoms=int(line)
        for i in range(num_atoms):
            line=file_unit.readline().split()
            temp_frame.coors.append(map(float,line[3:6]))
        line=file_unit.readline().split()
        temp_frame.box[0][0]=10.0*float(line[0])
        temp_frame.box[1][1]=10.0*float(line[1])
        temp_frame.box[2][2]=10.0*float(line[2])
        temp_frame.coors=array(temp_frame.coors,order='Fortran')
        temp_frame.coors=10.0*temp_frame.coors
        temp.append(temp_frame)
        pos=file_unit.tell()
        test_line=file_unit.readline()
        if len(test_line)==0: end=0

    return temp,0

def read_next(file_unit,iovars=None,iopos=None):

    temp_frame=cl_frame()
    file_unit.seek(pos)
    line=file_unit.readline()
    line=file_unit.readline()
    num_atoms=int(line)
    for i in range(num_atoms):
        line=file_unit.readline().split()
        temp_frame.coors.append(map(float,line[3:6]))
    line=file_unit.readline().split()
    temp_frame.box[0][0]=10.0*float(line[0])
    temp_frame.box[1][1]=10.0*float(line[1])
    temp_frame.box[2][2]=10.0*float(line[2])
    temp_frame.coors=array(temp_frame.coors,order='Fortran')
    temp_frame.coors=10.0*temp_frame.coors
    temp.append(temp_frame)
    pos=file_unit.tell()

    return temp_frame,pos,0

def read_frame(file_unit,iovars=None,iopos=None):

    return None,None,True

def close_traj(file_unit):

    file_unit.close()
    return None,0
    

def write_all():

    pass

def write_frame():

    pass
