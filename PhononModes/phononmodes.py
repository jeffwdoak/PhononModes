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
    - num_atom_types - number of atomic species in the crystal
    - atom_types - numpy array of the number of atoms of each type, shape:
          [num_atom_types]
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
        self.read_masses(mass_name)

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
        for i in range(self.num_modes):
            normal_modes.append([])
            in_file.readline()  # Discard blank line
            for j in range(self.num_atoms):
                normal_modes[i].append([])
                line = in_file.readline().split()
                for k in range(3):
                    normal_modes[i][j].append(float(line[k+1]))
        in_file.close()
        self.normal_modes = np.array(normal_modes)

    def read_masses(self,mass_name):
        """
        Reads in the atomic masses from the apos.dat file.
        """
        mass_file = open(mass_name,"r")
        # Discard unit cell parameters.
        for i in range(3):
            mass_file.readline()
        # Get number of unit cells
        line = mass_file.readline().split()
        num_cells = int(line[0])*int(line[1])*int(line[2])
        # Get number of atom types
        self.num_atom_types = int(mass_file.readline().split()[0])
        # Get list of atomic masses (strings).
        masses = mass_file.readline().split()
        # Get list containing number of atoms of each type (strings).
        atom_types = mass_file.readline().split()
        self.atom_types = np.array( [ int(i)*num_cells for i in atom_types ] )
        #self.atom_types = np.zeros(self.num_atom_types)
        #for i in atom_types:
        #    self.atom_types[i] = int(atom_types[i]) 
        mass_file.close()
        mass_vec = []
        for i in range(self.num_atom_types):
            for j in range(self.atom_types[i]):
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
        This definition is given by Eq. (1) of Pailhes, et al., PRL 113, 025506 
        (2014). This version corrects an error with the definition given in Eq.
        (4) of Hafner and Krajci, J. Phys.: Condens. Matter 5 (1993) 2489, where
        the term 'square' is not squared after the summation, leading to values
        of the participation ratio which are not normalized between 0 and 1.
        """
        square = 0.0
        quad = 0.0
        for j in range(len(self.normal_modes[index])):
            temp = np.dot(self.normal_modes[index,j],self.normal_modes[index,j])
            temp = temp/self.masses[j]
            square += temp
            quad += temp**2
        part = square**2/quad/self.num_atoms
        return part

    def atomic_participation_ratio(self,index):
        """
        Function to calculate the atom-decomposed participation ratio for the
        j-th atom in the k-th normal mode of a crystal. The atomic participation
        ratio is taken from Eq. (2) of Pailhes, et al., PRL 113, 025506 (2014).
        I do not believe that the atomic-participation-ratio is normalized
        between 0 and 1. Function returns a numpy array of length self.num_atoms
        containing the atomic participation ration of each atom in the mode
        """
        part = np.zeros(self.num_atoms)
        quad = 0.0
        for j in range(self.num_atoms):
            temp = np.dot(self.normal_modes[index,j],self.normal_modes[index,j])
            temp = temp/self.masses[j]
            part[j] = temp
            quad += temp**2
        quad = np.sqrt(self.num_atoms*quad)
        part = part/quad
        return part

    def list_of_atoms_by_species(self):
        """
        Function to sort each atom by its atomic species. Returns a nested list
        containing each atom, shape: [num_atom_types, atom_# ].
        E.g.,
        [[0, 1, 2, 3, 4, 5],
        [6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15]] for a 16-atom crystal, wherein atoms 1-6
        are of type A, atoms 7-8 are of type B, and atoms 9-16 are of type C.
        """
        species = [ [] for i in self.atom_types ]
        j = 0
        for i in range(self.num_atoms):
            if i >= np.sum(self.atom_types[:j+1]):
                j += 1
            print i,j
            species[j].append(i)
        return species


