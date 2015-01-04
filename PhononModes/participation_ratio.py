#!/usr/bin/env python

from phononmodes import *
import numpy as np
import sys


modes = PhononModes()

print "Mode_# Frequency_(cm^-1) Participation_Ratio"
for i in range(modes.num_modes):
    print i,modes.freqs[i],modes.participation_ratio(i)

