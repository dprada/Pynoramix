#from pynoramix_beta import *
from pynoramix import *
import pylab as pylab
import cPickle as cpickle
import mis_funcs as mi

#prot=molecule('aux.pdb',with_bonds=True)
prot2=molecule('1b23.pdb',with_bonds=True)

#hbonds_pack(system1=prot2,system2=prot2)
#aaaa=hbonds_pack(system1=prot2,system2=prot2)
#for ii in range(25):

#aa=hbonds_pack(system1=prot2)
for ii in range(50):
    hbonds(system1=prot2)

#aa=hbonds_pack(system1=prot2,system2=prot2)
#for ii in range(50):
#    hbonds(system1=prot2,system2=prot2,pack=aa)

#hbonds(system1=prot2,system2=prot2)
#hbonds_pack(system1=prot2,system2=prot2)
