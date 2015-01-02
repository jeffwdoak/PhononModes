#!/usr/bin/env python

from phononmodes import *
import numpy as np
import sys


modes = PhononModes()

index = int(sys.argv[1])

pr = modes.participation_ratio(index)

print pr

