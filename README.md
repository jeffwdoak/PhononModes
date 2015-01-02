PhononModes
===========

Class to read in and store phonon mode eigenvalues, eigenvectors, and atomic
masses from DMsymm-formatted phonons.out and apos.dat files.

Object attributes:
- num_atoms
- num_modes
- freqs - numpy array of phonon eigenvalues,shape: [num_modes]
- normal_modes - numpy array of phonon eigenvectors, shape: [num_modes, num_atoms, 3]
- masses - numpy array of atomic masses, shape: [num_atoms]

Object methods:
- participation_ratio:
    Function to calculate the participation ratio of the k-th normal mode of a
    crystal, where k is given by the argument index. The participation ratio is
    a value for a mode that ranges between 0 and 1, and gives a qualitative
    measure of the localization of a phonon mode: a value of 0 is localized in
    space while a value of 1 is delocalized across the entire crystal.
