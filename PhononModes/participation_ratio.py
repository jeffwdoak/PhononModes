#!/usr/bin/env python

"""
Script to calculate and output to screen the participation ratios of a phonon
calculation.

Assumes the phonon mode file is named 'phonons.out', the atomic mass file is
named 'apos.dat', and both are located in the directory where this script is
run.

"""

from phononmodes import *

modes = PhononModes()

print "Mode_# Frequency_(cm^-1) Participation_Ratio"
for i in range(modes.num_modes):
    print i, modes.freqs[i], modes.participation_ratio(i)
