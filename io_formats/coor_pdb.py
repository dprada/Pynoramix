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
    io_vars=[]
    io_pos=0
    io_err=False
    FFF=open(file_name,'r')
    return FFF,io_vars,io_pos,io_err # io_file,io_vars,io_pos,io_err

def read_all(file_unit,io_vars=None,io_pos=None):

    io_err=False
    io_end=False

    model_inds=[]
    file_unit.seek(0)
    while True:
        line=file_unit.readline()
        if line.startswith('MODEL'):
            model_inds.append([int(line.split()[1]),file_unit.tell()])
        if len(line)==0:
            break

    #pos=file_unit.tell()
    #file_unit.seek(pos)

    if len(model_inds)==0: model_inds.append([1,0])

    file_unit.seek(0)
    temp=[]
    for ref_mod in model_inds:
        frame=cl_frame()
        file_unit.seek(ref_mod[1])
        while True:
            line=file_unit.readline()
            ii=line.split()
            if ii[0]=='CRYST1':
                frame.box[0][0]=float(ii[1])
                frame.box[1][1]=float(ii[2])
                frame.box[2][2]=float(ii[3])
            if (ii[0] in ['ATOM','HETATM']):
                aux=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
                frame.coors.append(aux)
            if ii[0] in ['MODEL','END']:
                break

        frame.coors=array(frame.coors,order='Fortran')
        temp.append(frame)

    io_end=True

    return temp,io_err,io_end # io_file,io_err,io_end

def read_next(file_unit,iopos=None):

    return None,None,True,False  # io_file,io_pos,io_err,io_end

def read_frame(file_unit,iopos=None):

    return None,None,1

def close_traj(file_unit):

    io_err=False
    file_unit.close()
    return io_err  #io_err
    

def write_all():

    pass

def write_frame():

    pass
