####################################
# GENERAL COMMENTS

NAME_VERSION="Pynoramix 0.1"

#
# This module requires:
#

### External libraries:
from os import system
from os import path
from os import sys
import copy
import pylab
from numpy import *
import pickle as pic
import datetime as datetime

### pyno libraries:
from pyn_cl_coors import *
import top_par as tp
import pyn_fort_general as f
import pyn_math as pyn_math

#
# Structure of the file:
#
# -Class Set
#            - Instantiation
#            - Functions
# -Class Unit
# -Class Residue
# -Class Water
# -External functions
#

#
# END GENERAL COMMENTS
####################################


#######################################################
#######################################################
#### CLASSES
    
####
#### Labels or common attributes to be inherited
####

class labels_unit():                           # Every unit (atom) has at least these attributes
    def __init__(self):
        self.name=None
        self.index=None
        self.pdb_index=None
        self.covalent_bonds=[]                 # List of atoms covalently bonded.
        self.hbonds=[]                         # Atom h-bonded: [[atom_index][strength]]
        self.type=None                         # for residues: protein,ion,solv. for atom: atom_nature (H,O,...)

class labels_set(labels_unit):                 # Every set of units (chain, residue, molecule) has at least these attributes
    def __init__(self):
        self.num_atoms=0
        self.list_atoms=[]

class labels_parent(labels_unit):               # A selection always keep memory of the previous location
    def __init__(self,parent,argument=None,unit=False):
        if unit:
            self.name=parent.name
            self.index=parent.index
            self.pdb_index=parent.pdb_index
        else:
            self.name=parent.name
            self.condition=argument

####
#### Class Unit (atom)
####


class cl_unit(labels_unit):                     # Attributes of an atom

    def __init__(self):
        '''Initialize an atom object'''

        # From labels_unit: .name, .index, .pdb_index, .covalent_bonds    

        # > Topological properties

        self.resid=labels_unit()        # Residue which this atom belongs to. 
        self.chain=labels_unit()        # Chain which this atom belongs to.
        self.parent=None                # Selection which this atom comes from.

        self.alt_loc=0                  # Alternate location (PDB)
        self.code_ins_res=0             # Code of insertion of residues (PDB)
        self.seg_ident=''               # Index segment (PDB)
        self.elem_symb=''               # Element symbol (PDB)
        self.type_pdb=''                # Type of atom for the PDB (ATOM,HETATM)

        self.covalent_bonds=[]          # esto deberia estar heredado de labels_unit @!!!!!!@

        # > Physical and Chemical properties

        self.coors=[]
        self.mass=0.0                   # Mass
        self.charge=0.0                 # Charge
        self.vdw=0.0                    # VdW radius
        self.occup=0.0                  # Occupation (PDB)
        self.bfactor=0.0                # B-Factor
        self.acceptor=False             # True or false 
        self.donor=False                # True or false
        self.polar_class='Nothing'      # Acceptor or Donnor or nothing
        self.polarizability=False       # True of falsel

####
#### Class residue (set of atoms)
####


class cl_residue(labels_set):           # Attributes of a residue (inherits from label_set)

    def __init__( self ):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        self.chain=labels_unit()         # Chain which this residue belongs to.

        pass


####
#### Class water (set of atoms)
####

class cl_water(labels_set):             # Attributes of a water molecule
    '''Water specific attributes and methods'''

    def __init__( self ,water_model=None):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        self.model=water_model
                
        self.O=labels_unit()
        self.H1=labels_unit()
        self.H2=labels_unit()
        self.uvect_norm=[]
        self.microstate=''
        #def get_uvect_norm(self):
        
        pass

####
#### Class set (set of atoms: Molecule)
####


class molecule(labels_set):               # The suptra-estructure: System (waters+cofactors+proteins...)

    
    def __init__(self,input_file=None,download=None,coors=True,verbose=True,with_bonds=False):

        # From labels_set: .name, .index, .pdb_index, .num_atoms, .list_atoms

        # > Instantation options:
        self.file_topol=input_file       # pdb,gro...
        self.file_hbonds=''             # still not useful -do not remove-
        self.file_mss=''                # still not useful -do not remove-
        self.file_shell=''              # still not useful -do not remove-
        #self.selection=False            # is the set a selection (a mirror -pointer- of a parent set)
        
        # > Topological properties
        self.atom=[]                    # list of atoms    (objects: cl_unit)
        self.resid=[]                   # list of residues (objects: molecule)
        self.chain=[]                   # list of chains   (objects: molecule)
        self.chains=[]                  # list of chain names (strings)
        self.ion=[]                     # list of ions (objects: molecule)
        self.water=[]                   # list of waters   (objects: cl_water)
        self.water_model=None           # water model
        self.parent=None                # Parent set if it comes from selection (object: labels_parent)

        # > Physical and Chemical properties
        self.acceptors=[]               # list of acceptors
        self.donors=[]                  # list of donnors
        self.donors_hydrogen={}         # dict. (index donor atom : index H)
        self.dimensionality=0           # dimensionality (num_atoms*3)

        # > Coordinates and Trajectory
        self.which_traj=None
        self.which_frame=None
        self.box=[]
        self.cell=[]
        self.traj=[]
        self.dist_matrix=[]             # distance matrix (This should go to cl_coors?)

        # > Info PDB
        self.pdb_header=[]              # PDB Header (HEAD + TITLE)
        self.pdb_ss=[]                  # PDB Secondary structure
        

        ##################################

        # A SET CAN BE BUILT FROM A FILE OR FROM A SELECTION

        # IF IT COMES FROM A FILE, DOES IT NEED TO BE DOWNLOADED?

        if download:
            if not download.endswith('.pdb'):
                    download=download+'.pdb'
            input_file=download
            if not path.exists(input_file):
                temp='wget -nv http://www.rcsb.org/pdb/files/'+input_file+' 1>/dev/null 2>&1'
                system(temp)
                if path.exists(input_file):
                    print '# File saved as '+input_file
                else:
                    print '# Error downloading http://www.rcsb.org/pdb/files/'+input_file

            else:
                print '# The file '+input_file+' exists in the local folder. Loading it...'

            self.file_topol=input_file

        # BUILDING THE SET FROM A FILE

        if input_file:

            # The file does not exist?
            if not download and not path.exists(input_file):
                print "# The file "+input_file+" does not exist."
                return

            # Reading the file and
            # attaching the atoms: self.atom[] (cl_unit)

            self.load_topol(self.file_topol)

            # Finnal set up of the attributes of the cl_units in self.atom[]:
            
            before_resid=None
            before_chain=None
            jj=-1
            ii=-1
            kk=-1
            for atom in self.atom:
                if atom.resid.pdb_index!=before_resid :
                    before_resid=atom.resid.pdb_index
                    jj+=1   
                if atom.chain.name!=before_chain :
                    before_chain=atom.chain.name
                    kk+=1
                ii+=1
                atom.index=ii                   # atom index for pynoramix
                atom.resid.index=jj             # resid index for pynoramix
                atom.chain.index=kk
                atom.hbonds=[]
                atom.type=tp.atom_nature[tp.atom[atom.name]]
                atom.resid.type=tp.residue_type[atom.resid.name]
            ### Setting up the subsets: residues, waters, chains.

            # Auxiliary dictionary to build resids and waters.

            aux={}

            for atom in self.atom[:]:
                ii=atom.resid.index
                try: 
                    aux[ii][atom.name]=atom.index
                except:
                    aux[ii]={}
                    aux[ii][atom.name]=atom.index

            # Residues:

            for ii in aux.keys():
                temp_residue=cl_residue()
                temp_residue.index=ii
                temp_residue.list_atoms=aux[ii].values()
                jj=temp_residue.list_atoms[0]
                temp_residue.pdb_index=self.atom[jj].resid.pdb_index
                temp_residue.name=self.atom[jj].resid.name
                temp_residue.type=tp.residue_type[temp_residue.name]
                self.resid.append(temp_residue)


            # Waters

            for residue in self.resid[:]:
                if tp.residue_type[residue.name]=='Water':
                    if self.water_model==None:
                        self.water_model='tip'+str(len(residue.list_atoms))+'p'
                    temp_water=cl_water()
                    temp_water.name=residue.name
                    temp_water.index=residue.index
                    temp_water.pdb_index=residue.pdb_index
                    temp_water.list_atoms=residue.list_atoms
                    temp_water.model=self.water_model

                    for ii in residue.list_atoms:
                        atom=self.atom[ii]
                        if tp.atom[atom.name] in ['atOW']:
                            atom.polar_class='acceptor'
                            atom.polarizability=True
                            temp_water.O.index=atom.index
                            temp_water.O.pdb_index=atom.pdb_index
                            temp_water.O.name=atom.name
                        if tp.atom[atom.name] in ['atHW1','atHW2']:
                            atom.polar_class='donor'
                            atom.polarizability=True
                            if tp.atom[atom.name] in ['atHW1']:
                                temp_water.H1.index=atom.index
                                temp_water.H1.pdb_index=atom.pdb_index
                                temp_water.H1.name=atom.name
                            elif tp.atom[atom.name] in ['atHW2']:
                                temp_water.H2.index=atom.index
                                temp_water.H2.pdb_index=atom.pdb_index
                                temp_water.H2.name=atom.name

                    self.water.append(temp_water)

            # Chains
            ii=-1
            for atom in self.atom[:]:
                if atom.type_pdb in ['ATOM']:
                    if atom.chain.name not in self.chains:
                        ii+=1
                        self.chains.append(atom.chain.name)
                        temp_chain=labels_set()
                        temp_chain.name=atom.chain.name
                        temp_chain.index=ii
                        self.chain.append(temp_chain)
                    self.chain[ii].list_atoms.append(atom.index)
            for residue in self.resid[:]:
                ii=residue.list_atoms[0]
                residue.chain.name=self.atom[ii].chain.name
                residue.chain.index=self.atom[ii].chain.index

            # Ions

            for atom in self.atom[:]:
                if tp.residue_type[atom.resid.name]=='Ion':
                    temp_residue=cl_residue()
                    temp_residue.list_atoms=[atom.index]
                    temp_residue.index=atom.resid.index
                    temp_residue.pdb_index=atom.resid.pdb_index
                    temp_residue.name=atom.resid.name
                    self.ion.append(temp_residue)

            # Deleting the auxiliary dictionary:
            del(aux)

            ### Setting up the local attributes

            # Topology and Covalent bonds

            missing_atoms=[]
            if with_bonds:
                for residue in self.resid[:]:
                    # Missing atoms in the residue topology
                    found_atoms={}
                    for ii in residue.list_atoms:
                        found_atoms[tp.atom[self.atom[ii].name]]=ii
                    for ii in tp.residue_atoms[residue.name]:
                        if ii not in found_atoms.keys():
                            missing_atoms.append([ii,residue.name,residue.pdb_index])
                    # Covalent bonds
                    #for ii in tp.covalent_bonds[residue.name]:
                    #    aa=aux_name[ii[0]]
                    #    bb=aux_name[ii[1]]
                    #    self.atom[aa].covalent_bonds.append(bb)
                    #    self.atom[bb].covalent_bonds.append(aa)

            # Charge

            for atom in self.atom[:]:
                if tp.atom[atom.name] in tp.charge:
                    atom.charge=tp.charge[tp.atom[atom.name]]

            # Acceptors-Donors

                # Default:
            for atom in self.atom[:]:
                if tp.atom[atom.name] in tp.donors: atom.donor=True
                if tp.atom[atom.name] in tp.acceptors: atom.acceptor=True

                # Exceptions: (This needs to be polished)
            exc_res_don=[ii[0] for ii in tp.donors_exception]
            exc_res_acc=[ii[0] for ii in tp.acceptors_exception]
            for residue in self.resid[:]:
                if residue.name in exc_res_don:
                    for exception in tp.donors_exception:
                        if residue.name == exception[0]:
                            for ii in residue.list_atoms:
                                if tp.atom[self.atom[ii].name]==exception[1]:
                                    cov=0
                                    for jj in self.atom[ii].covalent_bonds:
                                        if tp.atom_nature[self.atom[jj].name]=='H': 
                                            cov=1
                                            break
                                        if exception[2]=='Always': self.atom[ii].donor=exception[2]
                                        elif exception[2]=='Hbonded' and cov==1: self.atom[ii].donor=exception[2]
                                        elif exception[2]=='Not Hbonded' and cov==0: self.atom[ii].donor=exception[2]
                if residue.name in exc_res_acc:
                    for exception in tp.acceptors_exception:
                        if residue.name == exception[0]:
                            for ii in residue.list_atoms:
                                if tp.atom[self.atom[ii].name]==exception[1]:
                                    cov=0
                                    for jj in self.atom[ii].covalent_bonds:
                                        if tp.atom_nature[self.atom[jj].name]=='H': 
                                            cov=1
                                            break
                                        if exception[2]=='Always': self.atom[ii].acceptor=exception[2]
                                        elif exception[2]=='Hbonded' and cov==1: self.atom[ii].acceptor=exception[2]
                                        elif exception[2]=='Not Hbonded' and cov==0: self.atom[ii].acceptor=exception[2]

            for atom in self.atom[:]:
                if atom.acceptor: self.acceptors.append(atom.index)
                if atom.donor:
                    self.donors.append(atom.index)
                    for ind_cov in atom.covalent_bonds:
                        if tp.atom_nature[tp.atom[self.atom[ind_cov].name]]=='H':
                            try: 
                                self.donors_hydrogen[atom.index].append(ind_cov)
                            except:
                                self.donors_hydrogen[atom.index]=[ind_cov]

            self.acceptors=array(self.acceptors,order='Fortran')
            self.donors=array(self.donors,order='Fortran')

            ### Setting up the global attributes

            self.name=self.file_topol[:-self.file_topol[::-1].find('.')-1]       # file=FOO.N.pdb -> name=FOO.N
            self.num_atoms=len(self.atom)
            self.dimensionality=self.num_atoms*3
            self.num_residues=len(self.resid)
            self.num_waters=len(self.water)
            self.num_chains=len(self.chains)
            self.num_ions=len(self.ion)
            self.list_atoms=[ii for ii in range(self.num_atoms)]

            ### Loading coordinates
            if coors:
                self.load_traj(self.file_topol,frame='ALL',verbose=False)


            if verbose:
                for ii in missing_atoms:
                    print '#',ii[0],ii[1],ii[2]
                self.info()
                if coors:
                    self.traj[0].info(index=0)

        ## END of IF input_file


    # END OF INSTANTATION

    ###
    ### INTERNAL FUNCTIONS OF A SET
    ###

    # Info function

    def info(self):

        self.num_atoms=len(self.atom)
        self.num_residues=len(self.resid)
        self.num_chains=len(self.chain)
        self.num_waters=len(self.water)
        self.num_trajs=len(self.traj)
        self.num_ions=len(self.ion)
        print '#','System created from the file',self.file_topol,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',self.num_chains,' chains'
        print '#',self.num_waters,' waters'
        print '#',self.num_ions,' ions'

    # To handle files

    def load_topol (self,name_file):

        if name_file.endswith('pdb'):
            ff=open(name_file,'r')
            for line in ff:
                ss=line.split()
                if ss[0] in ['HEADER','TITLE','CRYST1']: self.pdb_header.append(line)
                if ss[0].startswith('END'): break  # To read only the 1st model
                if ss[0] in ['HELIX','SHEET','TURN']: self.pdb_ss.append(line)
                if ss[0] in ['ATOM','HETATM']:
                    temp_atom=cl_unit()
                    temp_atom.type_pdb=line[0:6].replace(' ', '')
                    temp_atom.pdb_index=int(line[6:11])
                    temp_atom.name=(line[12:16].split())[0]
                    temp_atom.alt_loc=line[16]
                    temp_atom.resid.name=(line[17:20]).replace(' ', '')
                    temp_atom.chain.name=line[21]
                    temp_atom.resid.pdb_index=int(line[22:26])
                    temp_atom.code_ins_res=line[26]
                    temp_atom.occup=float(line[54:60])
                    temp_atom.bfactor=float(line[60:66])
                    temp_atom.seg_ident=line[72:76].replace(' ', '')
                    temp_atom.elem_symb=line[76:78].replace(' ', '')
                    temp_atom.charge=line[78:80].replace(' ', '')
                    temp_atom.index=len(self.atom)
                    self.atom.append(temp_atom)

            ff.close()

        if name_file.endswith('gro'):

            ff=open(name_file,'r')
            line=ff.readline()                                          # Header of the gro file
            line=ff.readline()                                        
            self.num_atoms=int(line)

            for i in range(self.num_atoms):           
            ## Fixed format taken from http://manual.gromacs.org/online/gro.html
                temp_atom=cl_unit()
                line=ff.readline()
                temp_atom.pdb_index=int(line[15:20])
                temp_atom.name=line[10:15].replace(" ", "")
                temp_atom.resid.name=line[5:10].replace(" ", "")
                temp_atom.resid.pdb_index=int(line[0:5]) 
                temp_atom.index=i           
                self.atom.append(temp_atom)

            ff.close()

    #def write_pdb (self,filename=None):
    #    
    #    if filename==None:
    #        print 'Enter filename: '
    #        print '      foo.write_pdb("foo.pdb")'
    #    else:
    #        if not filename.endswith('.pdb'): filename+='.pdb'
    #        if path.exists(filename): 
    #            print '# The file '+filename+' already exists.'
    #            return
    # 
    #        file=open(filename,'w')
    # 
    #        a='HEADER    '+'> CREATED BY PYNORAMIX '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+' <\n'
    #        file.write(str(a))
    # 
    #        for ii in self.pdb_header: file.write(str(ii))
    #        
    #        for ii in self.pdb_ss:     file.write(str(ii))
    # 
    #        dct_aux={'ATOM': 'ATOM  ', 'HETATM': 'HETATM'}
    #        
    #        new_index=0
    #        for ii in range(self.num_atoms):
    #            new_index+=1
    #            a=dct_aux[self.atom[ii].type_pdb]              # 1-6
    #            a+="%5d" % (new_index)                         # 7-11
    #            #a+="%5d" % self.atom[ii].pdb_index            # 7-11
    #            a+=' '                                         # 12
    #            a+=' '+"%-3s" % self.atom[ii].name             # 13-16
    #            a+=' '                                         # 17
    #            a+="%3s" % self.atom[ii].resid.name            # 18-20
    #            a+=' '                                         # 21
    #            a+="%1s" % self.atom[ii].chain.name            # 22
    #            a+="%4d" % self.atom[ii].resid.pdb_index       # 23-26
    #            a+=' '                                         # 27
    #            a+='   '                                       # 28-30
    #            a+="%8.3f" % float(self.frame[0].coors[ii][0])   # 31-38
    #            a+="%8.3f" % float(self.frame[0].coors[ii][1])   # 39-46
    #            a+="%8.3f" % float(self.frame[0].coors[ii][2])   # 47-54
    #            a+="%6.2f" % self.atom[ii].occup               # 55-60
    #            a+="%6.2f" % self.atom[ii].bfactor             # 61-66
    #            a+='          '                                # 67-76
    #            a+="%2s" % self.atom[ii].elem_symb             # 77-78
    #            a+="%2s" % self.atom[ii].charge                # 79-80
    #            a+='\n' 
    #            file.write(str(a))         
    #            if ii<(self.num_atoms-1):
    #                if self.atom[ii].type_pdb!=self.atom[ii+1].type_pdb or self.atom[ii].chain.name!=self.atom[ii+1].chain.name :
    #                    new_index+=1
    #                    a="TER   "
    #                    a+="%5d" % (new_index)
    #                    a+=' '
    #                    a+='  '
    #                    a+=' '                                         
    #                    a+="%3s" % self.atom[ii].resid.name            
    #                    a+='\n' 
    #                    file.write(str(a))
    #        a='END   '+'\n'
    #        file.write(str(a))
    #        file.close()
    #    return None

    def write_set_to_file(self,name_of_file):
        file=open(name_of_file,'w')
        pic.dump(self,file)
        file.close()

    def read_set_from_file(self,name_of_file):
        file=open(name_of_file,'r')
        A=pic.load(file)
        file.close()
        return A


    # To handle coordinates

    def info_trajs(self):
        if len(self.traj):
            for aa in range(len(self.traj)):
                self.traj[aa].info(index=aa)
        else:
            print '# No coordinates'
        pass

    def coors2frame (self,traj=0,frame=0):

        self.which_traj=traj
        self.which_frame=frame
        self.box=self.traj[traj].frame[frame].box
        for ii in range(self.num_atoms):
            self.atom[ii].coors=self.traj[traj].frame[frame].coors[ii,:]

        pass

    def load_traj (self,file_input=None,frame=None,begin=None,end=None,increment=1,units=None,verbose=True):

        temp_traj=cl_traj(file_input,frame,begin,end,increment,units,verbose=False)
        if verbose:
            temp_traj.info(index=len(self.traj))
        self.traj.append(temp_traj)
        del(temp_traj)
        if len(self.traj)==1:
            if len(self.traj[0].frame) :
                self.coors2frame()

    def delete_traj (self,index='ALL'):

        if index in ['all','All','ALL']:
            for ii in self.traj:
                if ii.io_opened:
                    ii.close()
            del(self.traj);self.traj=[]
            del(self.box);self.box=[]
            for ii in range(self.num_atoms):
                del(self.atom[ii].coors); self.atom[ii].coors=[]

        elif type(index) in ['int']:
            self.traj.__delitem__(index)
            if self.which_traj==index:
                del(self.box);self.box=[]
                for ii in range(self.num_atoms):
                    del(self.atom[ii].coors)
                    self.atom[ii].coors=[]
            elif self.which_traj>index:
                self.which_traj-=1

        pass

    def selection (self,condition=None,traj=0,frame='ALL',pbc=True):

        list_condition=selection(self,condition,traj,frame,pbc)
        return list_condition


###############################################################
###############################################################
    # To handle the set

    def distance(self,setA='ALL',setB=None,traj=0,frame='ALL',pbc=True):
        
        diff_system=1
        if pbc:
            check_cell=self.traj[traj].frame[0].cell
            if check_cell[0,1]!=90 or check_cell[0,2]!=90 or check_cell[1,2]!=90:
                print '# PBC not implemented for not orthorhombic boxes'
                return
            pbc=1

        if setA=='ALL' and setB in [None,'ALL','All','all']:
            diff_system=0

        if setA in ['ALL','All','all']:
            setA=[ii in range(self.num_atoms)]
            n_A=self.num_atoms
            natoms_A=self.num_atoms
        else:
            n_A=len(setA)
            natoms_A=self.num_atoms

        if setB in [None,'ALL','All','all']:
            setB=[ii in range(self.num_atoms)]
            n_B=self.num_atoms
            natoms_B=self.num_atoms
        else:
            n_B=len(setB)
            natoms_B=self.num_atoms

        if frame in ['ALL','All','all']:
            ll=len(self.traj[traj].frame)
            dists=zeros(shape=(n_A,n_B,ll),order='Fortran')
            for ii in range(ll):
                frame=self.traj[traj].frame[ii]
                dists[:,:,ii]=f.dist(diff_system,pbc,setA,frame.coors,frame.box,setB,frame.coors,n_A,n_B,natoms_A,natoms_B)


        if ll==1:
            return dists[:,:][0]
        else:
            return dists
            
    def rdf(self,setA=None,setB=None,traj=0,frame='ALL',pbc=True,bins=100,segment=None):

        if setA==None or setB==None:
            return

        if type(setA) not in [int32,int,list,tuple]:
            print type(setA)
            setA=self.selection(setA)
        elif type(setA) in [int,int32]:
            setA=[setA]

        if type(setB) not in [int32,int,list,tuple]:
            setB=self.selection(setB)
        elif type(setB) in [int,int32]:
            setB=[setB]

        n_A=len(setA)
        n_B=len(setB)
        natoms_A=self.num_atoms
        natoms_B=self.num_atoms

        if frame==0:
            frame=self.traj[traj].frame[frame]
            dist_frame=f.dist(1,pbc,setA,frame.coors,frame.box,setB,frame.coors,n_A,n_B,natoms_A,natoms_B)
            rdf_frame=f.rdf_frame(dist_frame,frame.box,segment[0],segment[1],bins,n_A,n_B)
            return rdf_frame

        elif frame in ['ALL','All','all']:
            xxx=pyn_math.binning(None,bins,segment,None,None)
            rdf_tot=zeros(shape=(bins),dtype=float,order='Fortran')
            num_frames=0.0
            for frame in self.traj[traj].frame:
                dist_frame=f.dist(1,pbc,setA,frame.coors,frame.box,setB,frame.coors,n_A,n_B,natoms_A,natoms_B)
                rdf_frame=f.rdf_frame(dist_frame,frame.box,segment[0],segment[1],bins,n_A,n_B)
                rdf_tot+=rdf_frame
                num_frames+=1.0
            rdf_tot=rdf_tot/num_frames
            return xxx,rdf_tot

        

    def neighbs(self,setA=None,setB=None,ranking=1,dist=None,traj=0,frame=0,pbc=True):
     
        if setA==None or setB==None:
            return

        if type(setA) not in [int,list,tuple]:
            setA=self.selection(setA)
        elif type(setA) in [int]:
            setA=[setA]

        if type(setB) not in [int,list,tuple]:
            setB=self.selection(setB)
        elif type(setB) in [int]:
            setB=[setB]

        n_A=len(setA)
        n_B=len(setB)
        natoms_A=self.num_atoms
        natoms_B=self.num_atoms

        diff_system=1

        if dist==None:
            
            if frame==0:
                frame=self.traj[traj].frame[frame]
                neighbs=f.neighbs_ranking(diff_system,pbc,ranking,setA,frame.coors,frame.box,setB,frame.coors,n_A,n_B,natoms_A,natoms_B)
                return neighbs
        
        else:
            if frame==0:
                frame=self.traj[traj].frame[frame]
                contact_map,num_neighbs,dist_matrix=f.neighbs_dist(diff_system,pbc,dist,setA,frame.coors,frame.box,setB,frame.coors,n_A,n_B,natoms_A,natoms_B)
                neighbs=[]
                if ranking:
                    for ii in range(n_A):
                        if num_neighbs[ii]:
                            neighbs_A=f.translate_list(1,setB,contact_map[ii,:],dist_matrix[ii,:],num_neighbs[ii],n_B)
                            neighbs.append(neighbs_A)
                        else:
                            neighbs.append([])
                else:
                    for ii in range(n_A):
                        if num_neighbs[ii]:
                            neighbs_A=f.translate_list(0,setB,contact_map[ii,:],dist_matrix[ii,:],num_neighbs[ii],n_B)
                            neighbs.append(neighbs_A)
                        else:
                            neighbs.append([])
                return neighbs

        print 'Not Implemented' 
        pass
        
    #def plot_contact_map(contact_map):
    #    
    #    pylab.gray()
    #    pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
    #    #pylab.matshow(contact_map==False)
    #    return pylab.show()
    # 
    #def rms_fit(self,set_reference=None,selection='all',new=False):
    #    
    #    coors_original=make_selection(self,selection).frame[0].coors
    #    coors_reference=make_selection(set_reference,selection).frame[0].coors
    # 
    #    if len(coors_original)!=len(coors_reference):
    #        print '# Error: Different number of atoms'
    #        return
    #    
    #    rot,center_ref,center_orig,rmsd,g=f.aux_funcs_general.min_rmsd(coors_reference,coors_original,len(coors_original))
    #    
    #    coors_original=self.frame[0].coors
    #    coors_new=f.aux_funcs_general.rot_trans(coors_original,rot,center_orig,center_ref,len(coors_original))
    # 
    #    if new:
    #        fitted_set=copy.deepcopy(self)
    #        fitted_set.frame[0].coors=coors_new
#   #         fitted_set.pdb_header="Mensaje del fitteo"
    #        fitted_set.rmsd=rmsd
    # 
    #        return fitted_set
    #    else:
    #        self.frame[0].coors=copy.deepcopy(coors_new)
    #        self.rmsd=rmsd
    # 
    #    print '# RMSD fitting:',rmsd
    #        # Use coors_new
    # 
    #    return
    # 
    #def displ_vector(self,set_reference=None):
    # 
    #    self.d_vector=set_reference.frame[0].coors - self.frame[0].coors





#### END CLASSES
#######################################################
#######################################################



#######################################################
#######################################################
#### EXTERNAL OBJECTS AND FUNCTIONS
    
####
#### Functions
####

#def min_distance(system,set_a,set_b=None,pbc=True,type_b='atoms'):
# 
#    if set_b==None:
# 
#        ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,True,system.frame[0].coors,system.frame[0].box,set_a,set_a,system.num_atoms,len(set_a),len(set_a))
# 
#    else:
#        if type_b=='atoms':
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms(pbc,False,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),len(set_b))
#        elif type_b=='vectors':
#            l_vects=shape(set_b)
#            if len(l_vects)==1:
#                set_b=[set_b]
#                l_vects=shape(set_b)
#            array(set_b,order='Fortran')
#            ind_a1,ind_a2,min_dist = f.aux_funcs_general.min_dist_atoms_ref(pbc,system.frame[0].coors,system.frame[0].box,set_a,set_b,system.num_atoms,len(set_a),l_vects[0])
#    
#    return ind_a1,ind_a2,min_dist
# 
# 
#def xtc2bin(xtc_name,bin_name):
#    command=home_path+'xtc2bin %s %s'%(xtc_name,bin_name)
# 
#    if path.exists(bin_name):
#        command2='mv %s %s#'%(bin_name,bin_name)
#        print 'file',bin_name,' was moved to ', bin_name+'#'
#        system(command2)
# 
#    system(command)
# 
# 
#def dot_product_3d(vect1,vect2):
# 
#    return f.aux_funcs_general.proj3d(vect1,vect2,len(vect1))
# 
#def isothermal_compressibility(system,Temp,input_file,frame=None,begin=None,end=None):
# 
#    frame=[ii for ii in range(begin,end+1)]
# 
#    V2a=0.0
#    Va=0.0
#    Kt=0.0
#    for ii in frame :
#        system.delete_coors()
#        system.load_coors (input_file,frame=ii)
#        xx=0.0
#        xx=system.frame[0].box[0][0]*system.frame[0].box[1][1]*system.frame[0].box[2][2]
#        Va+=xx
#        V2a+=xx**2
#    
#    Va=Va/(len(frame)*1.0)
#    V2a=V2a/(len(frame)*1.0)
#    Kt=(V2a-Va**2)/(Temp*Va)
#    
#    return Kt*0.10,'(nm/Kb)'



######################################################


def selection(system=None,condition=None,traj=0,frame='ALL',pbc=True):

    icondition=condition

    # attributes syntaxis:

    dict_selects={
        'backbone': '(atom.name N CA C O)',
        'sidechain': '(atom.resid.type Protein and not atom.name N CA C O H1 H2)',
        }

    for ii,jj in dict_selects.iteritems():
        icondition=icondition.replace(ii,jj)

    ### First block.
    icondition=icondition.replace(',',' ')
    icondition=icondition.replace('(',' ( ')
    icondition=icondition.replace(')',' ) ')
    icondition=icondition.replace('[',' ( ')
    icondition=icondition.replace(']',' ) ')
    icondition=icondition.replace('AND','and')
    icondition=icondition.replace('And','and')
    icondition=icondition.replace('OR','or')
    icondition=icondition.replace('Or','or')
    icondition=icondition.replace('NOT','not')
    icondition=icondition.replace('Not','not')
    icondition=icondition.replace('Within','within')
    icondition=icondition.replace('WITHIN','within')
    icondition=icondition.replace('Of','of')
    icondition=icondition.replace('OF','of')
    icondition=icondition.replace(' in ',' ')
    icondition=icondition.replace(' In ',' ')
    icondition=icondition.replace(' IN ',' ')
    ### Second block.
    icondition=icondition.replace(' chain',' atom.chain')
    icondition=icondition.replace(' resid',' atom.resid')

    # logic syntaxis

    icondition=icondition.split()
    aux_cond=[True for ii in range(len(icondition))]
    ocondition=[]

    for ii in range(len(icondition)):
        if icondition[ii] in ['(',')','and','or','not','within','of']:
            aux_cond[ii]=False
        if icondition[ii].startswith('atom'):
            aux_cond[ii]=False
    aux_cond.append(False)

    for ii in range(len(icondition)):
        part=icondition[ii]
        if part in ['(',')','and','or','not','of']:
            ocondition.append(part)
        elif part in ['within']:
            part2=icondition[ii+1]
            aux_cond[ii+1]=False
            ocondition.append(part)
            ocondition.append(float(part2))
        elif part.startswith('atom'):
            ocondition.append(part)
            if aux_cond[ii+1]:
                ocondition.append('in'); ocondition.append('[')
                jj=ii
                while True:
                    jj+=1
                    if not aux_cond[jj]:
                        ocondition.append(']')
                        break
                    else:
                        part2=icondition[jj]
                        try:
                            kk=float(part2)
                            ocondition.append(part2+',')
                        except:
                            ocondition.append("'"+part2+"',")
                        aux_cond[jj]=False
        else:
            if aux_cond[ii]:
                ocondition.append(part)
                aux_cond[ii]=False

    # Solving the 'withins':
    for ii in range(len(ocondition)):
        if 'within' == ocondition[ii]:
            cutoff=float(ocondition[ii+1])
            if ')' not in ocondition[0:ii]:
                sel1=[]
                cond1=' '.join(ocondition[0:ii])
                for atom in system.atom:
                    if eval(cond1):
                        sel1.append(atom.index)
            if '(' not in ocondition[(ii+3):]:
                sel2=[]
                cond2=' '.join(ocondition[(ii+3):])
                for atom in system.atom:
                    if eval(cond2):
                        sel2.append(atom.index)
            dists_sels=system.distance(sel1,sel2,traj=traj,frame=frame,pbc=pbc)
            list_sel1=[]
            for gg in range(len(dists_sels)):
                list_sel_frame=[]
                for jj in range(len(sel1)):
                    if f.within(dists_sels[gg][jj,:],cutoff,len(sel2)):
                        list_sel_frame.append(sel1[jj])
                list_sel1.append(list_sel_frame)
            #'(atom.resid.type Water and atom.type O) within 3.0 of atom.resid.type Protein'
            #atom.name OW within 3.0 of [3,4,5]
            #atom.name OW within 3.0 of sel1
            #atom.name OW within 3.0 of atom.name HW1
            if len(list_sel1)==1:
                return list_sel1[0]
            else:
                return list_sel1

    # Applying selection
    list_select=[]
    condition=' '.join(ocondition)    

    for atom in system.atom:
        if eval(condition):
            list_select.append(atom.index)
     
     
    return list_select










