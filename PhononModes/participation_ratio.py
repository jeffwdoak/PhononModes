#!/usr/bin/env python

from phononmodes import *
import numpy as np
import sys


modes = PhononModes()

print "Number of atoms:",modes.num_atoms
print "Number of modes:",modes.num_modes
print
print "Mode_# Frequency_(cm^-1) Participation_Ratio"
for i in range(modes.num_modes):
    print i,modes.freqs[i],modes.participation_ratio(i)

