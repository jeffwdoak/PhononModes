#!/usr/bin/env python

from phononmodes import *
import numpy as np
import sys


modes = PhononModes()

index = int(sys.argv[1])

aprs = modes.atomic_participation_ratio(index)

print "Number of atoms:",modes.num_atoms
print "Number of modes:",modes.num_modes
print "Participation Ratio from Eq. (1):"
print modes.participation_ratio(index)
print "Participation Ratio from Sum over Eq. (2):"
print aprs.sum()**2
print
print "Mode_# Atom_# Atomic_Participation_Ratio"
for j in range(modes.num_atoms):
    print index,j+1,aprs[j]

