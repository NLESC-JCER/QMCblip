from ase.calculators.calculator import FileIOCalculator, CalculatorSetupError
from ase.calculators.calculator import InputError
import numpy as np
from ase.units import Ang, Bohr, Ha
from os import path

class CHAMP(FileIOCalculator):

    implemented_properties = ['energy', 'forces', 'stress']

    name = "CHAMP"

    default_parameters = dict(
        vmc_in='vmc.inp',
        vmc_out='vmc.out',
        force_file='champ_forces',
        pos_file='molecule.xyz',
        champ_loc='/usr/bin/vmc.mov1',
        nodefile='',
        ncore=1)

    command = ""

    def __init__(self, restart=None,
                ignore_bad_restart_file=FileIOCalculator._deprecated,
                label='CHAMP', atoms=None, **kwargs):
        
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                label, atoms, **kwargs)

        if (not path.exists(self.parameters['champ_loc'])):
            raise CalculatorSetupError("Did not find champ at: " + self.parameters['champ_loc'])
            
        if (not path.exists(self.parameters['vmc_in'])):
            raise CalculatorSetupError("Problem reading CHAMP input.")

        if (self.parameters['nodefile'] == ""):
            self.command = "mpirun -n " + str(self.parameters['ncore']) + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " > " + self.parameters['vmc_out']
        else:
            self.command = "mpirun -s all -n " + str(self.parameters['ncore']) \
                    + " -machinefile " + self.parameters['nodefile'] + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " > " + self.parameters['vmc_out']


    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

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
            self.results['energy'] = float(fileobj.readline()) * Ha

            forces = []
            for line in fileobj.readlines():
                f_list = [float(i) for i in line.split(" ") if i.strip()]
                forces.append(f_list)

            self.results['forces'] = -np.asarray(forces) * Ha/Bohr   
        self.results['stress'] = np.zeros((6,)) 
