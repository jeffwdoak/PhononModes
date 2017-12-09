#!/usr/bin/env python

"""
Set of functions to compare phonon normal modes and frequencies between two
related structures to determine mode overlaps and anharmonic coupling constants.
These functions were used in Doak, Wolverton, Ozolins, Phys. Rev. B, 92, 174306
(2015) to calcuate Eq. 13--15.

When running this script, two command-line arguments are required:
    1. the path to a directory containing a phonon supercell calculation of the
    system at a saddle-point of a double-well mode.
    2. the path to a directory containing a phonon supercell calculation of the
    system at the bottom of double-well mode.

"""

import sys

import numpy as np

from phononmodes import *

saddledir = str(sys.argv[1])
welldir = str(sys.argv[2])

saddle = PhononModes(saddledir+"/phonons.out", saddledir+"/apos.dat")
well = PhononModes(welldir+"/phonons.out", welldir+"/apos.dat")

def test(i):
    """
    Calculate and compare participation ratios of modes with the same index in
    two calculations. Prints comparison to screen.

    Parameters
    ----------
    i : int
        Index of mode for which to calculate participation ratios in both
        calculations.

    """
    for j in range(well.num_modes):
        sp = saddle.participation_ratio(i)
        wp = well.participation_ratio(j)
        prod = np.dot(saddle.normal_modes[i].flatten(), well.normal_modes[j].flatten())
        print i, j, sp, wp, sp-wp, prod

def transformation_matrix(saddle, well):
    """
    Transformation matrix between two sets of normal mode eigenvectors.

    Write each normal mode of structure `saddle` as a linear combination of
    normal modes of the structure `saddle`: saddle_i = sum_j(A_ij*well_j)
    Transformation matrix, A_ij, gives the coefficients of this expansion.

    Parameters
    ----------
    saddle, well : PhononMode objects
        PhononModes for the saddle-point calculation and well-bottom
        calculations, respectively.

    Returns
    -------
    A : numpy array
        Transformation matrix between modes of saddle and modes of well.

    """
    A = np.dot(well.flattened_modes(), saddle.flattened_modes().T)
    return A

def D1(omega, omega_prime, sigma, initial_modes, final_modes):
    """
    Map frequencies from one structure onto another structure.

    From Vinay Hegde (hegdevinayi).
    Function from Vinay to map frequencies from one structure onto another
    structure.

    Parameters
    ----------
    omega, omega_prime : float
        Phonon frequencies to compare (units of cm**-1)
    sigma : float
        Tolerance below which to consider two frequencies as being similar.
    initial_modes, final_modes : PhononMode objects
        PhononModes for the two structures being compared.

    Returns
    -------
    D : float
        Overlap between modes with frequencies omega in structure 1 and modes
        with frequencies omega_prime in structure 2 (units of Angstroms).

    """
    D = 0.0
    for i in range(initial_modes.num_modes):
        for j in range(final_modes.num_modes):
            if np.abs(omega - initial_modes.freqs[i]) < sigma:
                if np.abs(omega_prime - final_modes.freqs[j]) < sigma:
                    D += (np.dot(initial_modes.normal_modes[i].flatten(),
                                 final_modes.normal_modes[j].flatten()
                                )**2)
    return D

def D0(mu, nu, initial, final):
    """
    Overlap between two modes in different structures.

    Calculates the overlap between mode index mu of structure intial and mode
    index nu of structure final.

    Parameters
    ----------
    mu, nu : int
        Indices of phonon modes in structures `intial` and `final`,
        respectively.
    initial, final : PhononMode objects
        PhononModes for the two structures to be compared.

    Returns
    -------
    D : float
        Overlap between modes (units of Angstroms).

    """
    D = np.dot(initial.normal_modes[mu].flatten(), final.normal_modes[nu].flatten())**2
    return D

def anharmonic_coupling0(well, saddle):
    """
    Anharmonic coupling constants between modes of two structures.

    Coupling constants are calculated from (Eq. 15 of Doak, PRB, 174306, 2015):
    a = (w_w**2 - w_s**2)*2/Pmax**2, Pmax=1

    Parameters
    ----------
    well, saddle : PhononMode objects
        PhononModes for the two structures to be compared.

    Returns
    -------
    a : numpy array
        Anharmonic coupling constants between modes of `well` and `saddle`
        (units of cm**-2). len(a) == saddle.num_modes

    """
    def square_freqs(freqs):
        square = np.zeros_like(freqs)
        for i in range(len(freqs)):
            if np.abs(freqs[i]) < 1.0:
                square[i] = 0.0
            elif freqs[i] < 0.0:
                square[i] = -(freqs[i]**2)
            else:
                square[i] = freqs[i]**2
        return square

    well_freqs_squared = square_freqs(well.freqs)
    saddle_freqs_squared = square_freqs(saddle.freqs)
    a = np.zeros(saddle.num_modes)
    for i in range(saddle.num_modes):
        a[i] = (well_freqs_squared[i] - saddle_freqs_squared[i])*2
    return a

def anharmonic_coupling1(well, saddle):
    """
    Anharmonic coupling constants between modes of two structures with acoustic
    modes removed.

    Coupling constants are calculated from (Eq. 15 of Doak, PRB, 174306, 2015):
    a = (w_w**2 - w_s**2)*2/Pmax**2, Pmax=1
    Coupling constants are not calculated for the acoustic modes of the two
    structures.

    Parameters
    ----------
    well, saddle : PhononMode objects
        PhononModes for the two structures to be compared.

    Returns
    -------
    a : numpy array
        Anharmonic coupling constants between modes of `well` and `saddle`
        (units of cm**-2). len(a) == saddle.num_modes

    """
    def square_freqs(freqs):
        square = np.zeros(len(freqs)-3)
        j = 0
        for i in range(len(freqs)):
            if np.abs(freqs[i]) < 1.0:
                continue
            elif freqs[i] < 0.0:
                square[j] = -(freqs[i]**2)
                j += 1
            else:
                square[j] = freqs[i]**2
                j += 1
        return square
    well_freqs_squared = square_freqs(well.freqs)
    saddle_freqs_squared = square_freqs(saddle.freqs)
    a = np.zeros(saddle.num_modes-3)
    for i in range(saddle.num_modes-3):
        a[i] = (well_freqs_squared[i] - saddle_freqs_squared[i])*2
    return a

def whittle_overlap(saddle, well, iter_):
    """
    Remove a number of modes with the highest overlap between two structures.

    Find and remove the saddle-point and well-bottom modes with the
    highest overlap (D0) from the overlap matrix, `iter_` times.

    Parameters
    ----------
    saddle, well : PhononMode objects
        PhononModes of the two structures to compare.
    iter_ : int
        Number of modes to identify and remove based on having the highest
        overlap between structuers.

    Returns
    -------
    overlap : numpy array
        Overlap between remaining modes of the two structures (units of
        Angstroms).
    welltosaddle : numpy array
        Array indicating which modes of `saddle` overlap most with modes of
        `well`. For each mode in `saddle`, this array contains the index of the
        mode in `well` that has the greatest overlap with it, up to `iter_`
        modes. If the mode in `saddle` does not have an overlap in the `iter_`
        highest overlaps, the value of its index in `welltosaddle` is 0.
        len(welltosaddle) == saddle.num_modes

    """
    # initial overlap matrix construction
    nmodes = saddle.num_modes
    saddlefreqs = np.copy(saddle.freqs)
    wellfreqs = np.copy(well.freqs)
    overlap = [D0(i, j, saddle, well) for i in range(nmodes) for j in range(nmodes)]
    overlap = np.reshape(overlap, (nmodes, nmodes))
    argoverlap = [[i, j] for i in range(nmodes) for j in range(nmodes)]
    argoverlap = np.reshape(argoverlap, (nmodes, nmodes, 2))
    welltosaddle = np.zeros(nmodes)
    # repeatedly find the maximum overlap and remove those two modes from the
    # saddle and well lists
    for i in range(iter_):
        max = np.max(overlap)
        index = np.unravel_index(np.argmax(overlap), np.shape(overlap))
        welltosaddle[argoverlap[index[0], index[1], 0]] = int(argoverlap[index[0], index[1], 1])
        print i, max, index, saddlefreqs[index[0]], wellfreqs[index[1]]
        overlap = np.delete(overlap, index[0], 0)
        overlap = np.delete(overlap, index[1], 1)
        argoverlap = np.delete(argoverlap, index[0], 0)
        argoverlap = np.delete(argoverlap, index[1], 1)
        saddlefreqs = np.delete(saddlefreqs, index[0])
        wellfreqs = np.delete(wellfreqs, index[1])

    return overlap, welltosaddle

def square_freqs(freqs):
    """
    Square phonon frequencies, accounting for imaginary modes.

    Squares each phonon frequency. Frequencies of imaginary modes are set to be
    negative. Square frequencies of acoustic modes are set to zero.

    Parameters
    ----------
    freqs : numpy array
        Array of phonon frequencies in cm**-1. Imaginary modes have negative
        frequencies.

    Returns
    -------
    square : numpy array
        Array of squared frequencies. Imaginary modes have negative square
        frequencies.
        len(square) == len(freqs)

    """
    square = np.zeros_like(freqs)
    for i in range(len(freqs)):
        if np.abs(freqs[i]) < 1.0:
            square[i] = 0.0
        elif freqs[i] < 0.0:
            square[i] = -(freqs[i]**2)
        else:
            square[i] = freqs[i]**2
    return square

def rotated_freqs(well, saddle):
    """
    Squared frequencies of modes of `well` acting on system `saddle`.

    Calculate the squared frequencies of the saddle system normal mode
    displacements acting on the dynamical matrix of the well system. Equation
    14 of Doak, PRB 174306, 2015.

    Parameters
    ----------
    well, saddle : PhononMode objects
        PhononModes of two structures to compare. Dynamical matrix of `well`
        will be used for comparison.

    Returns
    -------
    newfreqs : numpy array
        Array of effective square phonon frequencies (units of cm**-2) of modes
        of structure `saddle` acting on system `well`.

    """
    A = transformation_matrix(saddle, well)
    omegasquare = square_freqs(saddle.freqs)
    newfreqs = np.zeros_like(well.freqs)
    for i in range(len(well.freqs)):
        for j in range(len(saddle.freqs)):
            newfreqs[i] += omegasquare[j]*A[i, j]**2
        newfreqs[i] = np.sqrt(newfreqs[i])
    return newfreqs
