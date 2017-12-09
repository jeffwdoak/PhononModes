#!/usr/bin/env python

"""
Class to work with phonon normal modes calculated from the home-brew code
DMsymm. The class relies on two input files to populate its data structures.
These files are 'phonons.out' and 'apos.dat' The format of these files are shown
below.

phonons.out
-----------
This file contains the phonon frequencies (eigenvalues) and normal modes
(eigenvectors) calculated by DMsymm. The frequencies are given in cm**-1.

The first line of the file contains the # of atoms in the supercell and the # of
q-points in th calculation.
Then, for each q-point, there are the following sets of lines:
The next line contains the coordinate of the next q point in reciprocal space and
its weight.
The next 3n/6 lines contain the phonon frequencies of the q point, where n is
the number of atoms in the supercell. There are 6 frequencies per line.
The next line is blank.
The next 3*n*(n+1) lines contain normal mode displacements for each atom in each
normal mode, followed by single blank lines. A line in the mode displacement
section contains 
    (i) the index of the atom in the supercell,
    (ii) the real displacement along x, y, and z of the atom in the mode, and
    (iii) the imaginary displacement along x, y, and z of the atom in the mode.
There are n lines for each mode followed by a blank line. There are 3n modes
total.
After the last blank line, the file repeats this section with the next q point.

apos.dat
--------
This file contains information on the atomic crystal for which the phonons have
been calculated, the q-points at which the phonons have been calculated, and the
positions and masses of atoms in the crystal.

- Lines 1-3 - unit cell vectors of the crystal
- Line 4 - supercell size along x, y, and z
- Line 5 - # of atom types in the crystal and # of q-points in the calculation
- Line 6 - mass of each atom type
- Line 7 - number of atoms of each type in the unit cell (not the supercell)
- Line 8-n+8 - positions of each atom in cartesian coordinates of the unit cell
  (not the supercell)
- Line n+9 blank
- Line n+10-n+10+m - q-point coordinates in reciprocal space of the unit cell
  and q-point weights

Example apos.dat file:

0.0000 3.2814 3.2814
3.2814 0.0000 3.2814
3.2814 3.2814 0.0000
3 3 3
2 1
207.2 127.6
1 1
0.0000 0.0000 0.0000 Pb
3.2814 3.2814 3.2814 Te

0.0 0.0 0.0 1.0

"""

import sys

import numpy as np

from unitcell import *

class PhononModes(object):
    """
    Class to read in and store phonon mode eigenvalues, eigenvectors, and atomic
    masses from DMsymm-formatted phonons.out and apos.dat files.

    Attributes
    ----------
    num_atoms : int
        Number of atoms in the supercell.
    num_modes : int
        Number of phonon modes in the calculation.
    freqs : numpy array
        Array of phonon-mode frequencies (eigenvalues) in cm**-1. 
        len(freqs) == num_modes
    normal_modes : numpy array
        Array of phonon normal modes (eigenvectors).
        np.shape(normal_modes) == (num_modes, num_atoms, 3)
    num_atom_types : int
        Number of atomic species in the crystal.
    atom_types : numpy array
        Array of the number of atoms of each type in the supercell.
        len(atom_types) == num_atom_types
    masses : numpy array
        Array of atomic masses in g/mol. 
        len(masses) == num_atoms

    """

    def __init__(self, phon_name="phonons.out", mass_name="apos.dat"):
        """
        Create a new PhononModes object.

        Parameters
        ----------
        phon_name : str
            Name of file containing DMsymm-formatted eigenvalues and
            eigenvectors. Defaults to 'phonons.out'.
        mass_name : str
            Name of file contaiing atomic masses. Defaults to 'apos.dat'.

        """
        self.read_normal_modes(phon_name)
        self.read_masses(mass_name)

    def read_normal_modes(self, in_name):
        """
        Read in phonon normal mode coordinates from a phonons.out file.

        Parameters
        ----------
        phon_name : str
            Name of file containing DMsymm-formatted eigenvalues and
            eigenvectors. Defaults to 'phonons.out'.

        """
        with open(in_name, 'r') as in_file:
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
       self.normal_modes = np.array(normal_modes)

    def read_masses(self, mass_name):
        """
        Reads in the atomic masses from the apos.dat file.

        Parameters
        ----------
        mass_name : str
            Name of file contaiing atomic masses. Defaults to 'apos.dat'.

        """
        with open(mass_name, 'r') as mass_file:
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
            self.atom_types = np.array([int(i)*num_cells for i in atom_types])
       mass_vec = []
        for i in range(self.num_atom_types):
            for j in range(self.atom_types[i]):
                #for k in range(3):
                mass_vec.append(float(masses[i]))
        self.masses = np.array(mass_vec)

    def flattened_modes(self):
        """
        Flatten the array of normal modes.

        Returns
        -------
        flattened_modes : numpy array
            Array of normal modes with each mode being a 3*num_atoms long
            vector.

        """
        return np.array([i.flatten() for i in self.normal_modes])

    def participation_ratio(self, index):
        """
        Participation ratio of the k-th normal mode of a crystal.

        Calculate the participation ratio of the k-th normal mode of a crystal,
        where k is given by the argument index. The participation ratio is a
        value for a mode that ranges between 0 and 1, and gives a qualitative
        measure of the localization of a phonon mode: a value of 0 is localized
        in space while a value of 1 is delocalized across the entire crystal. 
        This definition is given by Eq. (1) of Pailhes, et al., PRL 113, 025506
        (2014). This version corrects an error with the definition given in Eq.
        (4) of Hafner and Krajci, J. Phys.: Condens. Matter 5 (1993) 2489, where
        the term 'square' is not squared after the summation, leading to values
        of the participation ratio which are not normalized between 0 and 1.

        Parameters
        ----------
        index : int
            Index of the normal mode array for which to calculate participation
            ratios.

        Returns
        -------
        part : float
            Participation ratio of mode `index`.

        """
        square = 0.0
        quad = 0.0
        for j in range(len(self.normal_modes[index])):  # sum over atoms
            temp = np.dot(self.normal_modes[index, j], self.normal_modes[index, j])
            temp = temp/self.masses[j]
            square += temp
            quad += temp**2
        part = square**2/quad/self.num_atoms
        return part

    def atomic_participation_ratio(self, index):
        """
        Atom-decomposed participation ratio of the k-th normal model of a
        crystal.

        Calculate the atom-decomposed participation ratio for the j-th atom in
        the k-th normal mode of a crystal. The atomic participation ratio is
        taken from Eq. (2) of Pailhes, et al., PRL 113, 025506 (2014).
        I do not believe that the atomic-participation-ratio is normalized
        between 0 and 1.

        Parameters
        ----------
        index : int
            Index of the normal mode for which to calculate atom-decomposed
            participation ratios.

        Returns
        -------
        part : numpy array
            Atom-decomposed participation ratios for each atom in mode `index`.

        """
        part = np.zeros(self.num_atoms)
        quad = 0.0
        for j in range(self.num_atoms):
            temp = np.dot(self.normal_modes[index, j], self.normal_modes[index, j])
            temp = temp/self.masses[j]
            part[j] = temp
            quad += temp**2
        quad = np.sqrt(self.num_atoms*quad)
        part = part/quad
        return part

    def list_of_atoms_by_species(self):
        """
        Sort each atom by its atomic species.
        
        Returns a nested list containing each atom, shape:
        [num_atom_types, atom_# ]
        e.g.,

        [[0, 1, 2, 3, 4, 5],
        [6, 7],
        [8, 9, 10, 11, 12, 13, 14, 15]]

        for a 16-atom crystal, wherein atoms 1-6 are of type A, atoms 7-8 are of
        type B, and atoms 9-16 are of type C.

        Returns
        -------
        species : nested list
            List containing a list for each type of atom. Within the secondary
            lists are the indices of each atom with that atom type.

        """
        species = [[] for i in self.atom_types]
        j = 0
        for i in range(self.num_atoms):
            if i >= np.sum(self.atom_types[:j+1]):
                j += 1
            species[j].append(i)
        return species
