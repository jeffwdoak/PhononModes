PhononModes
===========

Module to work with phonon normal modes calculated from the home-brew code
DMsymm written by Prof. Vidvuds Ozolins. The class relies on two input files to
populate its data structures. These files are 'phonons.out' and 'apos.dat'. 
The format of these files are shown below.

This module has been used to calculate low-temperature ferroelectric phase 
transitions in the semiconductor alloy system PbS--PbTe in
Doak, Wolverton, Ozolins, Physical Review B, 92, 174306 (2015).

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

