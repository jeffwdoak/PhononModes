#!/usr/bin/env python

"""
Script to calculate atom-decomposed participation ratios of a phonon
calculation.

Assumes the phonon mode file is named 'phonons.out', the atomic mass file is
named 'apos.dat', and both are located in the directory where this script is
run.

Script takes 1 command line argument, the index of the phonon mode for which to
calculate the atom-decomposed participation ratio.

Prints atom-decomposed participation ratio to standard output.

"""

import sys

from phononmodes import *

modes = PhononModes()

index = int(sys.argv[1])

aprs = modes.atomic_participation_ratio(index)

print "Number of atoms:", modes.num_atoms
print "Number of modes:", modes.num_modes
print "Participation Ratio from Eq. (1):"
print modes.participation_ratio(index)
print "Participation Ratio from Sum over Eq. (2):"
print aprs.sum()**2
print
print "Mode_# Atom_# Atomic_Participation_Ratio"
for j in range(modes.num_atoms):
    print index, j+1, aprs[j]
