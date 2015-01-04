#!/usr/bin/env python

# Script to plot the atomic participation ratio of each atom in each mode,
# separated by atomic species. The result is a number of lists of (frequency, 
# atomic particpation ratio) tuples, one for each type of atom in the crystal.

from phononmodes import *
import numpy as np
import sys


modes = PhononModes()

species_list = modes.list_of_atoms_by_species()


apr_lists = [ [] for i in modes.atom_types ]

for k in range(modes.num_modes):
    apr = modes.atomic_participation_ratio(k)
    apr_tuple=np.array([ [ modes.freqs[k] for i in apr ], apr]).T
    for i in range(modes.num_atom_types):
        apr_lists[i].append(apr_tuple[species_list[i]])
for i in range(len(apr_lists)):
    apr_lists[i] = np.array(apr_lists[i]).reshape((-1,2))

for i in range(len(apr_lists)):
    print "Frequency_(cm^-1) Participation_Ratio"
    for j in range(len(apr_lists[i])):
        print apr_lists[i][j,0],apr_lists[i][j,1]
    print


sys.exit()
