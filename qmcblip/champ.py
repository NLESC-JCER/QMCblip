"""ASE CHAMPS Calculator."""
from os import path

import numpy as np
from ase.calculators.calculator import (CalculatorSetupError, FileIOCalculator,
                                        InputError)
from ase.units import Ang, Bohr, Ha

from qmcblip.champio import Settings


class CHAMP(FileIOCalculator):
    """ASE Calculator for CHAMP.

    An ASE calculator for the Quantum Monte Carlo software CHAMP.

    Args:
        vmc_in (:obj:`str`, optional): File to read CHAMP settings from.
        vmc_out (:obj:`str`, optional):The output file of CHAMP.
        force_file (:obj:`str`, optional): The file that CHAMP writes the forces and energies to.
        pos_file (:obj:`str`, optional): The file from which CHAMP will read to
            location of the atoms.
        champ_loc (:obj:`str`, optional): Location of the CHAMP executable.
        nodefile (:obj:`str`, optional): If set, the calculator will run on
            multiple nodes for CHAMP.
        ncore (:obj:`str`, optional): Amount of cores to run CHAMP on.
        settings (:obj:`Settings<qmcblip.champio.Settings>`, optional): Input settings for CHAMP
            (OVERRULES VMC.INP).
        use_opt_wf (:obj:`bool`, optional): Use the optimized WF from last step.
    """

    # Stress only for FLARE compatibility
    implemented_properties = ['energy', 'forces', 'stress']

    name = "CHAMP"

    # The parameters to run CHAMP
    default_parameters = dict(
        vmc_in='vmc.inp',
        vmc_out='vmc.out',
        force_file='write_forces',
        pos_file='molecule.xyz',
        champ_loc='/usr/bin/vmc.mov1',
        nodefile=None,
        ncore=1,
        settings=None,
        use_opt_wf=False)

    # Placeholder
    command = ""

    def __init__(self, restart=None,
                ignore_bad_restart_file=FileIOCalculator._deprecated,
                label='CHAMP', atoms=None, **kwargs):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                label, atoms, **kwargs)

        if not path.exists(self.parameters['champ_loc']):
            raise CalculatorSetupError("Did not find champ at: " + self.parameters['champ_loc'])

        if (not path.exists(self.parameters['vmc_in']) and self.parameters['settings'] is None):
            raise CalculatorSetupError("You did no supply any configuration parameters for CHAMP. \
                                       Either give a .inp file or give a configuration dictionary.")

        if (not isinstance(self.parameters['settings'], Settings) and
            self.parameters['settings'] is not None):
            raise CalculatorSetupError('You did not provide a proper settings object.')

        if self.parameters['settings'] is not None:
            self.parameters['settings'].write(filename=self.parameters['vmc_in'])
        else:
            self.parameters['settings'] = Settings.read(self.parameters['vmc_in'])

        # Set the command to call CHAMP
        self._set_command()

    def _set_command(self):
        # Add nodefile or not
        if self.parameters['nodefile'] is None:
            self.command = "mpirun -n " + str(self.parameters['ncore']) + " " \
                    + self.parameters['champ_loc'] + " -i vmc_temp.inp" \
                    + " -o " + self.parameters['vmc_out']
        else:
            self.command = "mpirun -s all -n " + str(self.parameters['ncore']) \
                    + " -machinefile " + self.parameters['nodefile'] + " " \
                    + self.parameters['champ_loc'] + " -i vmc_temp.inp" \
                    + " -o " + self.parameters['vmc_out']

    def configure(self, **kwargs):
        """(Re)configure the CHAMP calculator by setting the keyword arguments.

        The tags for the vmc.inp can also be set here.

        Args:
            **kwargs: keyword arguments.
        """
        self.set(**kwargs)

        # Set the command to call CHAMP. May have changed due to the parameters.
        self._set_command()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        # Write the coordinates of the atoms in xyz format
        with open(self.parameters['pos_file'], 'w', encoding='utf-8') as fileobj:
            comment = atoms.symbols
            natoms = len(atoms)
            fmt='%22.15f'

            fileobj.write('%d\n%s\n' % (natoms, comment))

            for label, (x, y, z) in zip(atoms.symbols, atoms.positions):
                fileobj.write('%-2s %s %s %s\n' % (label, fmt % (x * Ang/Bohr), \
                                fmt % (y * Ang/Bohr), fmt % (z * Ang/Bohr)))

        if self.parameters['use_opt_wf']:
            self.parameters['settings'].use_opt_wf()

        # Rewrite the input file
        self.parameters['settings'].write('vmc_temp.inp')



    def read(self, label):
        raise NotImplementedError

    def read_results(self):
        if not path.exists(self.parameters['force_file']):
            raise InputError("Could not detect the forces from CHAMP in " \
                                        + self.parameters['force_file'])

        with open(self.parameters['force_file'], 'r', encoding="utf-8") as fileobj:
            # Read the energy and convert to eV
            self.results['energy'] = float(fileobj.readline()) * Ha

            # Read the forces
            forces = []
            for line in fileobj.readlines():
                f_list = [float(i) for i in line.split(" ") if i.strip()]
                forces.append(f_list)

            # Convert to eV/A and use ASE force direction
            self.results['forces'] = -np.asarray(forces) * Ha/Bohr

        # Pass an null stress tensor for FLARE compatibility
        self.results['stress'] = np.zeros((6,))
