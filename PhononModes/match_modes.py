#!/usr/bin/env python

from phononmodes import *
import numpy as np
import sys

saddledir = str(sys.argv[1])
welldir = str(sys.argv[2])

saddle = PhononModes(saddledir+"/phonons.out",saddledir+"/apos.dat")
well = PhononModes(welldir+"/phonons.out",welldir+"/apos.dat")


def test(i):
    for j in range(well.num_modes):
        sp = saddle.participation_ratio(i)
        wp = well.participation_ratio(j)
        prod = np.dot(saddle.normal_modes[i].flatten(),well.normal_modes[j].flatten())
        print i,j,sp,wp,sp-wp,prod


def transformation_matrix(saddle,well):
    """
    Function to determine the transformation matrix between two sets of normal
    mode eigenvectors.
    """
    A = np.dot(well.flattened_modes(),saddle.flattened_modes().T)
    return A

def D1(omega,omega_prime,sigma,initial_modes,final_modes):
    """
    Function from Vinay to map frequencies from one structure onto another
    structure.
    """
    D = 0.0
    for i in range(initial_modes.num_modes):
        for j in range(final_modes.num_modes):
            if np.abs(omega - initial_modes.freqs[i]) < sigma:
                if np.abs(omega_prime - final_modes.freqs[j]) < sigma:
                    D += np.dot(initial_modes.normal_modes[i].flatten(),final_modes.normal_modes[j].flatten())**2
    return D
                    
def D0(mu,nu,initial,final):
    """
    Calculates the overlap between mode index mu of structure intial and mode
    index nu of structure final.
    """
    D = np.dot(initial.normal_modes[mu].flatten(),final.normal_modes[nu].flatten())**2
    return D

def anharmonic_coupling0(well,saddle):
    """
    Function to determine anharmonic coupling constants a, from w_well and
    w_saddle: a = (w_w**2 - w_s**2)*2/Pmax**2, Pmax=1
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

def anharmonic_coupling1(well,saddle):
    """
    Function to determine anharmonic coupling constants a, from w_well and
    w_saddle: a = (w_w**2 - w_s**2)*2/Pmax**2, Pmax=1. This function also
    removes the acoustic modes from the summation coupling constants.
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

def whittle_overlap(saddle,well,iter):
    """
    Function to find and remove the saddle-point and well-bottom modes with the
    highest overlap (D0) from the overlap matrix, iter times.
    """
    # initial overlap matrix construction
    nmodes = saddle.num_modes
    saddlefreqs = np.copy(saddle.freqs)
    wellfreqs = np.copy(well.freqs)
    overlap = [ D0(i,j,saddle,well) for i in range(nmodes) for j in range(nmodes) ]
    overlap = np.reshape(overlap,(nmodes,nmodes))
    argoverlap = [ [i,j] for i in range(nmodes) for j in range(nmodes) ]
    argoverlap = np.reshape(argoverlap,(nmodes,nmodes,2))
    welltosaddle = np.zeros(nmodes)
    # repeatedly find the maximum overlap and remove those two modes from the
    # saddle and well lists
    for i in range(iter):
        max = np.max(overlap)
        index = np.unravel_index(np.argmax(overlap),np.shape(overlap))
        welltosaddle[argoverlap[index[0],index[1],0]] = int(argoverlap[index[0],index[1],1])
        print i,max,index,saddlefreqs[index[0]],wellfreqs[index[1]]
        overlap = np.delete(overlap,index[0],0)
        overlap = np.delete(overlap,index[1],1)
        argoverlap = np.delete(argoverlap,index[0],0)
        argoverlap = np.delete(argoverlap,index[1],1)
        saddlefreqs = np.delete(saddlefreqs,index[0])
        wellfreqs = np.delete(wellfreqs,index[1])

    return overlap,welltosaddle


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

def rotated_freqs(well,saddle):
    """
    Function to calculate the squared frequencies of the well system normal mode
    displacements acting on the dynamical matrix of the saddle system.
    """
    A = transformation_matrix(saddle,well)
    omegasquare = square_freqs(saddle.freqs)
    newfreqs = np.zeros_like(well.freqs)
    for i in range(len(well.freqs)):
        for j in range(len(saddle.freqs)):
            newfreqs[i] += omegasquare[j]*A[i,j]**2
        newfreqs[i] = np.sqrt(newfreqs[i])
    return newfreqs
