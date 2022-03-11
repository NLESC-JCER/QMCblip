from ase.calculators.calculator import FileIOCalculator, CalculatorSetupError
from ase.calculators.calculator import InputError
import numpy as np
from ase.units import Ang, Bohr, Ha
from os import path

class CHAMP(FileIOCalculator):
    """
    An ASE calculator for the Quantum Monte Carlo software CHAMP
    """ 

    # Stress only for FLARE compatibility
    implemented_properties = ['energy', 'forces', 'stress']

    name = "CHAMP"

    # The parameters to run CHAMP
    default_parameters = dict(
        vmc_in='vmc.inp',
        vmc_out='vmc.out',
        force_file='champ_forces',
        pos_file='molecule.xyz',
        champ_loc='/usr/bin/vmc.mov1',
        nodefile='',
        ncore=1)

    # Placeholder
    command = ""

    def __init__(self, restart=None,
                ignore_bad_restart_file=FileIOCalculator._deprecated,
                label='CHAMP', atoms=None, **kwargs):
        """
        Keyword arguments:
        vmc_in -- The input file for CHAMP (default 'vmc.inp')
        vmc_out -- The output file of CHAMP (default 'vmc.out')
        force_file -- The file that CHAMP writes the forces and energies to (default 'champ_forces')
        pos_file -- The file from which CHAMP will read to location of the atoms (default 'molecule.xyz')
        champ_loc -- Location of the CHAMP executable (default '/usr/bin/vmc.mov1')
        nodefile -- If set, the calculator will run on multiple nodes for CHAMP
        ncore -- Amount of cores to run CHAMP on (default 1)
        """
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                label, atoms, **kwargs)

        if (not path.exists(self.parameters['champ_loc'])):
            raise CalculatorSetupError("Did not find champ at: " + self.parameters['champ_loc'])
            
        if (not path.exists(self.parameters['vmc_in'])):
            raise CalculatorSetupError("Problem reading CHAMP input.")

        # Add nodefile or not
        if (self.parameters['nodefile'] == ""):
            self.command = "mpirun -n " + str(self.parameters['ncore']) + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " -o " + self.parameters['vmc_out']
        else:
            self.command = "mpirun -s all -n " + str(self.parameters['ncore']) \
                    + " -machinefile " + self.parameters['nodefile'] + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " -o " + self.parameters['vmc_out']


    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        # Write the coordinates of the atoms in xyz format
        with open(self.parameters['pos_file'], 'w') as fileobj:
            comment = atoms.symbols
            natoms = len(atoms)
            fmt='%22.15f'

            fileobj.write('%d\n%s\n' % (natoms, comment))

            for s, (x, y, z) in zip(atoms.symbols, atoms.positions):
                fileobj.write('%-2s %s %s %s\n' % (s, fmt % (x * Ang/Bohr), \
                                fmt % (y * Ang/Bohr), fmt % (z * Ang/Bohr)))

    def read(self, label):
        raise NotImplementedError

    def read_results(self):
        if (not path.exists(self.parameters['force_file'])):
            raise InputError("Could not detect the forces from CHAMP in " \
                                        + self.parameters['force_file'])

        with open(self.parameters['force_file'], 'r') as fileobj:
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
