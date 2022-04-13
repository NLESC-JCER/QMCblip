"""Tools for simulations with sGDML."""
import os
from os.path import exists

import numpy as np
from ase.io import read

from sgdml import __version__
from sgdml.utils import io

def db_to_sgdml(db_file, dataset_file, name=None):
    """Convert an ASE DB to a sGDML dataset.

    Note:
        Licensed under MIT license.
        Copyright (c) 2018-2020 Stefan Chmiela

    Args:
        db_file (:obj:`Path`): database filename.
        dataset_file (:obj:`Path`): sGDML dataset file.
        name (:obj:`str`): name of the new dataset.
    """

    # MIT License
    #
    # Copyright (c) 2018-2020 Stefan Chmiela
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in all
    # copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    # SOFTWARE.

    # File checks
    if not exists(db_file):
        raise FileNotFoundError(db_file + " is not in this directory!")
    if exists(dataset_file):
        print(dataset_file + " already exists. Overwriting...")

    basename = os.path.splitext(os.path.basename(db_file))[0]

    mols = read(db_file, index=':')

    # filter incomplete outputs from trajectory
    mols = [mol for mol in mols if mol.get_calculator() is not None]

    R, z, E, F = None, None, None, None

    calc = mols[0].get_calculator()

    if 'forces' not in calc.results:
        raise ValueError("Forces are missing in the database file!")

    Z = np.array([mol.get_atomic_numbers() for mol in mols])
    if not (Z == Z[0]).all():
        raise ValueError('Order of atoms changes accross database!')

    R = np.array([mol.get_positions() for mol in mols])
    z = Z[0]

    E = np.array([mol.get_potential_energy() for mol in mols])
    F = np.array([mol.get_forces() for mol in mols])

    if name is None:
        name = basename

    # Base variables contained in every model file.
    base_vars = {
        'type': 'd',
        'code_version': __version__,
        'name': name,
        'theory': 'QMC',
        'R': R,
        'z': z,
        'F': F,
    }

    base_vars['F_min'], base_vars['F_max'] = np.min(F.ravel()), np.max(F.ravel())
    base_vars['F_mean'], base_vars['F_var'] = np.mean(F.ravel()), np.var(F.ravel())

    base_vars['r_unit'] = 'Ang'
    base_vars['e_unit'] = 'eV'

    if E is not None:
        base_vars['E'] = E
        base_vars['E_min'], base_vars['E_max'] = np.min(E), np.max(E)
        base_vars['E_mean'], base_vars['E_var'] = np.mean(E), np.var(E)

    base_vars['md5'] = io.dataset_md5(base_vars)
    np.savez_compressed(dataset_file, **base_vars)
