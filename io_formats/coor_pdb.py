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
    return FFF,0,0

def read_all(file_unit,iopos=None):

    model_inds=[]
    file_unit.seek(0)
    for line in file_unit:
        if line.startswith('MODEL'):
            pos=file_unit.tell()
            model_inds.append(int(line.split()[1]),pos)


    #pos=file_unit.tell()
    #file_unit.seek(pos)

    if len(model_inds)==0: model_inds.append([1,0])

    temp=[]
    for ref_mod in model_inds:
        frame=cl_frame()
        file_unit.seek(ref_mod[1])
        for line in file_unit:
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

    return temp,0

def read_next(file_unit,iopos=None):

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
