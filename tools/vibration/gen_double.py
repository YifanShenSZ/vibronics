'''
According to user input generate vibration*.in,
specifying the symmetry adapted vibrational basis

The non-totally symmetric irreduciles are truncated to single and double excitations

User should write their input in section 'User input' in this script
'''

from typing import List, Tuple
from pathlib import Path
import copy
import numpy

''' User input '''
# point group product table, the totally symmetric irreducible must be 1
product_table = numpy.array([
   [1, 2],
   [2, 1]
])

NModes_per_irred = (27, 24)

max_phonons = (
    (1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 3, 0, 2, 1, 1, 1, 0, 0, 0, 1, 0),
    (0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)
''' End of user input '''

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
    # Single and double special: quick return if there are > 2 excited non-totally symemtric modes
    count_asym = 0
    for mode in excited_modes:
        if mode[0] != 0:
            count_asym += 1
            if count_asym > 2: return
    # quick return if there are no such excitations
    for mode in excited_modes:
        if max_phonons[mode[0]][mode[1]] == 0: return
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
        # Finish when counter overflows
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
    max_excitation = 0
    possible_modes = []
    excitation = 0
    for i in range(len(max_phonons)):
        for j in range(len(max_phonons[i])):
            if (max_phonons[i][j] > 0):
                max_excitation += 1
                possible_modes.append((i, j))
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
        # basic case: the leading `excitation` modes in `possible_modes` are excited
        excited_indices = [*range(excitation)]
        # Map indices to modes of irreducibles
        excited_modes = []
        for index in excited_indices:
            excited_modes.append(possible_modes[index])
        generate_all(excited_modes, files)
        # Loop over possible excited modes as an excitation-nary counter, with ascending digits
        while True:
            excited_indices[-1] += 1
            # Carry to former digit
            for i in range(-1, -excitation, -1):
                if excited_indices[i] >= max_excitation + 1 + i:
                    excited_indices[i - 1] += 1
                    excited_indices[i] = 0
            # Guarantee ascendance
            for i in range(excitation - 1):
                if excited_indices[i] >= excited_indices[i + 1]:
                    excited_indices[i + 1] = excited_indices[i] + 1
            # Finish when counter overflows
            if excited_indices[-1] >= max_excitation: break
            # Map indices to modes of irreducibles
            excited_modes = []
            for index in excited_indices:
                excited_modes.append(possible_modes[index])
            generate_all(excited_modes, files)
