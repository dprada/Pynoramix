from numpy import zeros
class cl_frame():
    def __init__(self):
        self.time=None
        self.step=None
        self.precision=None
        self.model=None
        self.coors=[]
        self.box=zeros(shape=(3,3),order='Fortran')
