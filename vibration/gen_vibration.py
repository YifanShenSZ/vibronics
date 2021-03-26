'''
According to user input generate vibration*.in,
specifying the symmetry adapted vibrational basis

User should write their input in section 'User input' in this script
'''

from typing import List, Tuple
from pathlib import Path
import copy
import numpy

''' User input '''
# point group product table, the totally symmetric irreducible must be 1
#product_table = numpy.array([
#   [1, 2],
#   [2, 1]
#])
#
#NModes_per_irred = (2, 1)
#
#max_phonons = (
#    (2, 2),
#    (2,)
#)

product_table = numpy.array([
   [1,],
])

NModes_per_irred = (3,)

max_phonons = (
    (9, 9, 9),
)
''' End of user input '''

# Given index of a mode
# Return its irreducible and its index in this irreducible
def index2mode(index:int) -> (int, int):
    index_copy = copy.deepcopy(index)
    for i in range(len(NModes_per_irred)):
        if index_copy < NModes_per_irred[i]: return (i, index_copy)
        index_copy -= NModes_per_irred[i]
    raise Exception("index out of number of normal modes")

def determine_excitation(phonons:List) -> int:
    excitation = 0
    for irred in phonons:
        for phonon in irred:
            if phonon > 0: excitation += 1
    return excitation

def determine_irreducible(phonons:List) -> int:
    irred = 0
    for i in range(1, len(phonons)):
        order = 0
        for phonon in phonons[i]:
            order += phonon
        if order % 2 == 1: irred = product_table[irred][i]
    return irred

# Given excited modes
# Return the basic case: the excited modes has |1>, others |0>
def generate_base(excited_modes:List) -> List:
    base = []
    for NModes in NModes_per_irred:
        base.append([])
        for i in range(NModes): base[-1].append(0)
    for mode in excited_modes:
        base[mode[0]][mode[1]] = 1
    return base

# Given excited modes
# Output all possible vibrations to files
def generate_all(excited_modes:List, files:List) -> List:
    # basic case
    phonons = generate_base(excited_modes)
    file = files[determine_irreducible(phonons)]
    for irred in phonons:
        for phonon in irred:
            print("%5d" % phonon, file=file, end='')
        print(file=file)
    # loop as a len(excited_modes)-nary counter
    while True:
        phonons[excited_modes[0][0]][excited_modes[0][1]] += 1
        # Carry to latter digit
        for i in range(1, len(excited_modes)):
            if phonons[excited_modes[i - 1][0]][excited_modes[i - 1][1]] \
            > max_phonons[excited_modes[i - 1][0]][excited_modes[i - 1][1]]:
                phonons[excited_modes[i - 1][0]][excited_modes[i - 1][1]] = 1
                phonons[excited_modes[i][0]][excited_modes[i][1]] += 1
        # break when counter overflows
        if phonons[excited_modes[-1][0]][excited_modes[-1][1]] \
        > max_phonons[excited_modes[-1][0]][excited_modes[-1][1]]: break
        file = files[determine_irreducible(phonons)]
        for irred in phonons:
            for phonon in irred:
                print("%5d" % phonon, file=file, end='')
            print(file=file)

if __name__ == "__main__":
    # Subtract 1 to product table for python convenience
    product_table -= 1
    # Sanity check
    NIrreds = len(NModes_per_irred)
    assert(NIrreds == len(max_phonons))
    for i in range(NIrreds): assert(NModes_per_irred[i] == len(max_phonons[i]))
    # Count
    NModes = sum(NModes_per_irred)
    max_excitation = determine_excitation(max_phonons)
    # Replace basis files
    files = []
    for i in range(NIrreds):
        files.append(open("vibration" + str(i + 1) + ".in", 'w'))

    # |0>
    ground_state = generate_base([])
    for irred in ground_state:
        for phonon in irred:
            print("%5d" % phonon, file=files[0], end='')
        print(file=files[0])

    # loop over excitation
    for excitation in range(1, max_excitation + 1):
        # Select basic excited modes by indices
        excited_indices = [*range(excitation)]
        # Map indices to modes of irreducibles
        excited_modes = []
        for index in excited_indices:
            excited_modes.append(index2mode(index))
        generate_all(excited_modes, files)
        # Loop over possible excited modes as an excitation-nary counter, with ascending digits
        while True:
            excited_indices[-1] += 1
            # Carry to former digit
            for i in range(-1, -excitation, -1):
                if excited_indices[i] >= NModes + 1 + i:
                    excited_indices[i - 1] += 1
                    excited_indices[i] = 0
            # Guarantee ascendance
            for i in range(excitation - 1):
                if excited_indices[i] >= excited_indices[i + 1]:
                    excited_indices[i + 1] = excited_indices[i] + 1
            # Finish when counter overflows
            if excited_indices[-1] >= NModes: break
            # Map indices to modes of irreducibles
            excited_modes = []
            for index in excited_indices:
                excited_modes.append(index2mode(index))
            generate_all(excited_modes, files)