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

    
    def __init__(self,input_file=None,download=None,coors=True,verbose=True,with_bonds=True,missing_atoms=True):

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
        self.donors=[[],[]]             # list of [donors,donor_hydrogen]
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

            # Residues:

            aux_dict={}
            for atom in self.atom[:]:
                ii=atom.resid.index
                try: 
                    aux_dict[ii][atom.name]=atom.index
                except:
                    aux_dict[ii]={}
                    aux_dict[ii][atom.name]=atom.index

            for ii in aux_dict.keys():
                temp_residue=cl_residue()
                temp_residue.index=ii
                temp_residue.list_atoms=aux_dict[ii].values()
                jj=temp_residue.list_atoms[0]
                temp_residue.pdb_index=self.atom[jj].resid.pdb_index
                temp_residue.name=self.atom[jj].resid.name
                temp_residue.type=tp.residue_type[temp_residue.name]
                self.resid.append(temp_residue)

            del(aux_dict)

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


            ### Setting up the local attributes

            # Topology and Covalent bonds

            aux_dict={}
            for atom in self.atom[:]:
                ii=atom.resid.index
                try: 
                    aux_dict[ii][tp.atom[atom.name]]=atom.index
                except:
                    aux_dict[ii]={}
                    aux_dict[ii][tp.atom[atom.name]]=atom.index

            if with_bonds:
                for residue in self.resid[:]:
                    jj=residue.index
                    tp_residue_name=tp.residue[residue.name]
                    tp_residue_atoms=tp.residue_atoms[tp_residue_name]
                    # Missing or unknown atoms in the residue topology
                    found_atoms=[]
                    for ii in residue.list_atoms:
                        found_atoms.append(tp.atom[self.atom[ii].name])
                    missing_atoms=set(tp_residue_atoms).difference(found_atoms)
                    unknown_atoms=set(found_atoms).difference(tp_residue_atoms)
                    # Covalent bonds: Topology residue
                    for ii in tp.covalent_bonds[tp_residue_name]:
                        try:
                            aa=aux_dict[jj][ii[0]]
                            bb=aux_dict[jj][ii[1]]
                            self.atom[aa].covalent_bonds.append(bb)
                            self.atom[bb].covalent_bonds.append(aa)
                        except:
                            pass
                    # Covalent bonds: Peptide bond [C(i)->N(i+1)]
                    try:
                        next_residue=self.resid[jj+1]
                        if residue.type=='Protein' and next_residue.type=='Protein':
                            if residue.chain.name==next_residue.chain.name :
                                aa=aux_dict[jj]['atC']
                                bb=aux_dict[jj+1]['atN']
                                self.atom[aa].covalent_bonds.append(bb)
                                self.atom[bb].covalent_bonds.append(aa)
                    except:
                        pass
                    # Covalent bonds: Terminals

                    # Listing missing or unknown atoms
                    if missing_atoms:
                        for ii in missing_atoms:
                            print '# No atom type',ii,'in', residue.name, residue.pdb_index
                        for ii in unknown_atoms:
                            print '# Unknown atom type',ii,'in', residue.name, residue.pdb_index

            del(aux_dict)

            # Charge

            for atom in self.atom[:]:
                if tp.atom[atom.name] in tp.charge:
                    atom.charge=tp.charge[tp.atom[atom.name]]

            # Acceptors-Donors

            for atom in self.atom[:]:
                # Donors default
                if tp.atom[atom.name] in tp.donors: 
                    atom.donor=True
                # Donors exceptions
                if tp.atom[atom.name] in tp.donors_with_exceptions:
                    try:
                        exception=donors_exception[atom.name][atom.resid.name]
                        if exception[0]=='Always':
                            atom.donor=exception[1]
                        elif exception[0]=='Without H':
                            if 'H' not in [self.atom[ii].type for ii in atom.covalent_bonds]:
                                atom.donor=exception[1]
                        else:
                            print '# Error with donors exceptions:', atom.name, atom.resid.name
                    except:
                        pass
                # Acceptors default
                if tp.atom[atom.name] in tp.acceptors: 
                    atom.acceptor=True
                # Donors exceptions
                if tp.atom[atom.name] in tp.acceptors_with_exceptions:
                    try:
                        exception=acceptors_exception[atom.name][atom.resid.name]
                        if exception[0]=='Always':
                            atom.acceptor=exception[1]
                        elif exception[0]=='With H':
                            if 'H' in [self.atom[ii].type for ii in atom.covalent_bonds]:
                                atom.acceptor=exception[1]
                        else:
                            print '# Error with acceptors exceptions:', atom.name, atom.resid.name
                    except:
                        pass
                
                if atom.acceptor:
                    self.acceptors.append(atom.index)

                if atom.donor:
                    for ii in atom.covalent_bonds:
                        if self.atom[ii].type=='H':
                            self.donors[0].append(atom.index)
                            self.donors[1].append(ii)

            self.acceptors=array(self.acceptors,dtype=int,order='Fortran')
            self.donors=array(self.donors,dtype=int,order='Fortran')

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
                    temp_atom.resid.name=(line[17:21]).replace(' ', '') # real format: 17:20
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

        elif type(index) in [int]:
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

    #def hbonds_selection(self,setA='ALL',setB=None,verbose=True):
    # 
    #    setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)
    #    num_donors_A=0
    #    num_accetp_A=0
    #    for ii in setA:
    #        if 
    #    num_donors_B=
    #    num_accept_B=

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

        setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)
        dists=empty(shape=(n_A,n_B,num_frames),order='Fortran')

        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            dists[:,:,num_frames]=f.dist(diff_system,pbc,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
            num_frames+=1

        if num_frames==1:
            return dists[:,:][0]
        else:
            return dists
            
    def rdf(self,setA=None,setB=None,traj=0,frame='ALL',pbc=True,bins=100,segment=None):

        setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)

        xxx=pyn_math.binning(None,bins,segment,None,None)
        rdf_tot=zeros(shape=(bins),dtype=float,order='Fortran')
        num_frames=0
        for iframe in __read_frame_opt__(self,traj,frame):
            dist_frame=f.dist(1,pbc,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
            rdf_frame=f.rdf_frame(dist_frame,frame.box,segment[0],segment[1],bins,n_A,n_B)
            rdf_tot+=rdf_frame
            num_frames+=1.0

        rdf_tot=rdf_tot/(num_frames*1.0)
        return xxx,rdf_tot


    def neighbs(self,setA=None,setB=None,ranking=1,dist=None,traj=0,frame=0,pbc=True):
     
        setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)
        num_frames=__length_frame_opt__(self,traj,frame)

        if dist==None:
            neighbs=[]
            neighbs=empty(shape=(n_A,ranking,num_frames),dtype=int,order='Fortran')
            num_frames=0
            for iframe in __read_frame_opt__(self,traj,frame):
                neighbs[:,:][num_frames]=f.neighbs_ranking(diff_system,pbc,ranking,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
                num_frames+=1
            if num_frames==1:
                return neighbs[:,:][0]
            else:
                return neighbs
            
        else:
            neighbs=[]
            sort_opt=0
            if ranking:
                sort_opt=1
            for iframe in __read_frame_opt__(self,traj,frame):
                contact_map,num_neighbs,dist_matrix=f.neighbs_dist(diff_system,pbc,dist,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
                aux_neighbs=[]
                if ranking:
                    for ii in range(n_A):
                        if num_neighbs[ii]:
                            neighbs_A=f.translate_list(sort_opt,setB,contact_map[ii,:],dist_matrix[ii,:],num_neighbs[ii],n_B)
                            aux_neighbs.append(neighbs_A)
                        else:
                            aux_neighbs.append([])
                neighbs.append(aux_neighbs)
            if num_frames==1:
                return neighbs[0]
            else:
                return neighbs

    #def contact_map (self,setA=None,setB=None,dist=None,traj=0,frame=0,pbc=True):
    # 
    #    setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,setA,None,setB)
    #    num_frames=__length_frame_opt__(self,traj,frame)
    # 
    #    contact_map=empty(shape=(n_A,n_B,num_frames),dtype=int,order='Fortran')
    #    num_frames=0
    #    for iframe in __read_frame_opt__(self,traj,frame):
    #        contact_map[:,:][num_frames]=f.contact_map(diff_system,pbc,dist,setA,iframe.coors,iframe.box,setB,iframe.coors,n_A,n_B,natoms_A,natoms_B)
    #        num_frames+=1
    # 
    #    if num_frames==1:
    #        return contact_map[:,:][0]
    #    else:
    #        return contact_map
    # 
    #    #def plot_contact_map(contact_map):
    #    #    
    #    #    pylab.gray()
    #    #    pylab.imshow(contact_map==False,origin='lower',interpolation=None) # Would be better not to interpolate
    #    #    #pylab.matshow(contact_map==False)
    #    #    return pylab.show()

    #def lrmsd_fit(self,set_ref=None,traj_ref=None,frame_ref=None,selection=None,traj=None,frame=None,new=False):
    # 
    #    setA,n_A,natoms_A,setB,n_B,natoms_B,diff_system=__read_sets_opt__(self,set_ref,None,selection)
    # 
    #    if n_A!=n_B :
    #        print '# Error: Different number of atoms'
    #        return
    # 
    #    rot,center_ref,center_orig,rmsd,g=f.aux_funcs_general.min_rmsd(coors_reference,coors_original,len(coors_original))
    #    coors_new=f.aux_funcs_general.rot_trans(coors_original,rot,center_orig,center_ref,len(coors_original))
    # 
    #    
    #    if new:
    #        fitted_set=copy.deepcopy(self)
    #        fitted_set.frame[0].coors=coors_new
    #         fitted_set.pdb_header="Mensaje del fitteo"
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

def __read_frame_opt__ (syst=None,traj=0,frame=None):

    #Options: integer,list,tuple,'ALL'

    if frame in ['ALL','All','all']:
        return syst.traj[traj].frame
    elif type(frame) in [int32,int]:
        return [syst.traj[0].frame[ii] for ii in [frame]]
    elif type(frame) in [list,tuple]:
        return [syst.traj[0].frame[ii] for ii in frame]

def __length_frame_opt__ (syst=None,traj=0,frame=None):

    #Options: integer,list,tuple,'ALL'

    if frame in ['ALL','All','all']:
        return syst.traj[traj].num_frames
    elif type(frame) in [int32,int]:
        return 1
    elif type(frame) in [list,tuple]:
        return len(frame)

def __read_sets_opt__(systA=None,setA=None,systB=None,setB=None):

    if setA==None:
        print '# SetA needed.'
        pass

    diff_system=1
    if systB==None:

        nsysa=systA.num_atoms
        nsysb=systA.num_atoms

        if setA in ['ALL','All','all']:
            setA=[ii in range(systA.num_atoms)]
            nlista=systA.num_atoms
        elif type(setA) in [int32,int]:
            setA=[setA]
            nlista=1
        elif type(setA) in [list,tuple]:
            nlista=len(setA)
        else:
            setA=systA.selection(setA)

        if setB in [None]:
            setB=[ii in range(systA.num_atoms)]
            nlistb=systA.num_atoms
        elif setB in ['ALL','All','all']:
            setB=[ii in range(systA.num_atoms)]
            nlistb=systA.num_atoms
        elif type(setB) in [int32,int]:
            setB=[setB]
            nlistb=1
        elif type(setB) in [list,tuple]:
            nlistb=len(setB)
        else:
            setB=systA.selection(setB)

        if (setA==setB): diff_system=0

        return setA,nlista,nsysa,setB,nlistb,nsysb,diff_system

    else:
        print 'Different systems: Not implemented yet'
        pass
          


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










