#!/usr/bin/env python

import numpy as np
from unitcell import *
import sys


class PhononModes(object):
    """
    Class to read in and store phonon mode eigenvalues, eigenvectors, and atomic
    masses from DMsymm-formatted phonons.out and apos.dat files.
    Object attributes:
    - num_atoms
    - num_modes
    - freqs - numpy array of phonon eigenvalues [num_modes]
    - normal_modes - numpy array of phonon eigenvectors, shape:
          [num_modes, num_atoms, 3]
    - masses - numpy array of atomic masses [num_atoms]
    Object methods:
    - participation_ratio
    """
    def __init__(self,phon_name="phonons.out",mass_name="apos.dat"):
        """
        Create a new phonon_modes object, with phonon eigenvalues and vectors
        read from file phon_name, and atomic masses read from file mass_name.
        """
        self.read_normal_modes(phon_name)
        self.read_atomic_masses(mass_name)

    def read_normal_modes(self,in_name):
        """
        Reads in phonon normal mode coordinates from a phonons.out file.
        """
        in_file = open(in_name,"r")
        line = in_file.readline().split()
        self.num_atoms = int(line[0])
        num_q_points = int(line[1])
        in_file.readline()  # Discard q-point line
        self.num_modes = 3*self.num_atoms
        num_lines = self.num_modes/6  # Six frequencies per line
        freqs = []
        for i in range(num_lines):
            line = in_file.readline().split()
            for j in line:
                freqs.append(float(j))
        self.freqs = np.array(freqs)
        normal_modes = []
        for i in range(num_modes):
            normal_modes.append([])
            in_file.readline()  # Discard blank line
            for j in range(num_atoms):
                normal_modes[i].append([])
                line = in_file.readline().split()
                for k in range(3):
                    normal_modes[i,j].append(float(line[k+1]))
        in_file.close()
        self.normal_modes = np.array(normal_modes)

    def read_masses(self,mass_name):
        """
        Reads in the atomic masses from the apos.dat file.
        """
        mass_file = open(mass_name,"r")
        # Discard unit cell parameters and supercell size.
        for i in range(4):
            mass_file.readline()
        # Get number of atom types
        num_atom_types = int(mass_file.readline().split()[0])
        # Get list of atomic masses (strings).
        masses = mass_file.readline().split()
        # Get list containing number of atoms of each type (strings).
        atom_types = mass_file.readline().split()
        mass_file.close()
        mass_vec = []
        for i in range(num_atom_types):
            for j in range(int(atom_types[i])):
                #for k in range(3):
                mass_vec.append(float(masses[i]))
        self.masses = np.array(mass_vec)

    def participation_ratio(self,index):
        """
        Function to calculate the participation ratio of the k-th normal mode of a
        crystal, where k is given by the argument index. The participation ratio is
        a value for a mode that ranges between 0 and 1, and gives a qualitative
        measure of the localization of a phonon mode: a value of 0 is localized in
        space while a value of 1 is delocalized across the entire crystal.
        """
        square = 0.0
        quad = 0.0
        for j in range(len(self.normal_modes[index])):
            temp = np.dot(self.normal_modes[index,j],self.normal_modes[index,j])
            square += temp/self.masses[j]
            quad += temp**2/self.masses[j]
        part = square/quad/self.num_atoms
        return part
