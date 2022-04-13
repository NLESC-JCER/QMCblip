"""Module for usefull tools."""
from os import remove
from os.path import exists

import numpy as np
import periodictable as pt
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.db import connect
from ase.io.trajectory import Trajectory

from flare import otf_parser

def traj_to_db(traj_file, db_file, append=False):
    """Convert ASE trajectory file to database.

    Additionally, the tag 'relaxed' is set to False.

    Args:
        traj_file (:obj:`Path`): trajectory filename.
        db_file (:obj:`Path`): database filename.
        append (:obj:`bool`): append to a current database if it exists.
            If database file does not exists this does nothing.
    """
    # File Checks
    if not exists(traj_file):
        raise FileNotFoundError(traj_file + " is not in this directory!")
    if exists(db_file) and append:
        print(db_file + " already exists. Appending to it...")
    elif exists(db_file) and not append:
        print(db_file + " already exists. Overwriting...")
        remove(db_file)

    # Open files
    database = connect(db_file)
    traj = Trajectory(traj_file)

    # Write to database
    for atoms in traj:
        database.write(atoms, relaxed=False)

def otf_to_db(otf_file, db_file, append=False):
    """Convert FLARE OTF file to database.

    Additionally, the tag 'relaxed' is set to False.

    Args:
        otf_file (:obj:`Path`): OTF filename.
        db_file (:obj:`Path`): database filename.
        append (:obj:`bool`): append to a current database if it exists.
            If database file does not exists this does nothing.
    """
    # File checks
    if not exists(otf_file):
        raise FileNotFoundError(otf_file + " is not in this directory!")
    if exists(db_file) and append:
        print(db_file + " already exists. Appending to it...")
    elif exists(db_file) and not append:
        print(db_file + " already exists. Overwriting...")
        remove(db_file)

    # Open files
    database = connect(db_file)
    otf = otf_parser.OtfAnalysis(otf_file, calculate_energy=True)

    types = otf.header['species']
    atomicnumbers = np.zeros(len(types))
    atomicmasses = np.zeros(len(types))
    for i, _ in enumerate(types):
        atomicnumbers[i] = getattr(pt, types[i]).number
        atomicmasses[i] = getattr(pt, types[i]).mass

    for i, _ in enumerate(otf.dft_frames):
        # Get all the interesting data
        forces = otf.gp_force_list[i]
        positions = otf.gp_position_list[i]
        energy = otf.gp_thermostat['potential energy'][i]

        # Set all the properties
        atoms = Atoms(''.join(types), positions=positions)
        atoms.calc = Calculator()
        atoms.calc.implemented_properties = ['forces', 'energy']
        atoms.calc.results['forces'] = forces
        atoms.calc.results['energy'] = energy
        atoms.numbers = atomicnumbers
        atoms.masses = atomicmasses

        database.write(atoms)
