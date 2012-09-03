from numpy import zeros
class cl_frame():
    def __init__(self):
        self.time=None
        self.step=None
        self.precision=None
        self.lam=None                    # Coming from trr files... what is this?
        self.model=None
        self.coors=[]
        self.box=zeros(shape=(3,3),order='Fortran')
        self.cell=zeros(shape=(3,3),order='Fortran')
