from ase.calculators.calculator import FileIOCalculator, CalculatorSetupError
from ase.calculators.calculator import InputError
import numpy as np
from ase.units import Ang, Bohr, Ha
from os import path, rename
from functools import reduce
from caf.champio import write_input, read_input, use_opt_wf

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
        force_file='write_forces',
        pos_file='molecule.xyz',
        champ_loc='/usr/bin/vmc.mov1',
        nodefile=None,
        ncore=1,
        tags=None,
        use_opt_wf=False)

    # Placeholder
    command = ""

    def __init__(self, restart=None,
                ignore_bad_restart_file=FileIOCalculator._deprecated,
                label='CHAMP', atoms=None, **kwargs):
        """
        Keyword arguments:
        vmc_in -- The input file for CHAMP (default 'vmc.inp')
        vmc_out -- The output file of CHAMP (default 'vmc.out')
        force_file -- The file that CHAMP writes the forces and energies to (default 'write_forces')
        pos_file -- The file from which CHAMP will read to location of the atoms (default 'molecule.xyz')
        champ_loc -- Location of the CHAMP executable (default '/usr/bin/vmc.mov1')
        nodefile -- If set, the calculator will run on multiple nodes for CHAMP
        ncore -- Amount of cores to run CHAMP on (default 1)
        tags -- Input tags for CHAMP
        use_opt_wf -- Use the optimized WF from last step (default False)
        """
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                label, atoms, **kwargs)

        if (not path.exists(self.parameters['champ_loc'])):
            raise CalculatorSetupError("Did not find champ at: " + self.parameters['champ_loc'])
            
        if (not path.exists(self.parameters['vmc_in'])):
            raise CalculatorSetupError("Problem reading CHAMP input.")

        if self.parameters['tags'] is not None:
            write_input(self.parameters['tags'], filename=self.parameters['vmc_in'])
        else:
            self.parameters['tags'] = read_input(self.parameters['vmc_in'])

        if self.parameters['use_opt_wf']:
            self.parameters['vmc_in'] = 'vmc_temp.inp'
            write_input(self.parameters['tags'], self.parameters['vmc_in'])

        self._set_command()
    
    def _set_command(self):
        # Add nodefile or not
        if (self.parameters['nodefile'] is None):
            self.command = "mpirun -n " + str(self.parameters['ncore']) + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " -o " + self.parameters['vmc_out']
        else:
            self.command = "mpirun -s all -n " + str(self.parameters['ncore']) \
                    + " -machinefile " + self.parameters['nodefile'] + " " \
                    + self.parameters['champ_loc'] + " -i " + self.parameters['vmc_in'] \
                    + " -o " + self.parameters['vmc_out']

    def configure(self, **kwargs):
        """
        (Re)configure the CHAMP calculator by setting the keyword arguments.
        The tags for the vmc.inp can also be set here.
        """
        self.set(**kwargs)

        if self.parameters['tags'] is not None:
            write_input(self.parameters['tags'], filename=self.parameters['vmc_in'])

        self._set_command()

    def configure_qmc(self, update_tags):
        if type(update_tags) is not dict:
            raise TypeError("You did not supply a dictionary to update the QMC settings!")

        tags = self.parameters['tags']
    	
        for key, item in update_tags.items():
            key = key.split('-')
            reduce(dict.__getitem__, key[:-1], tags)[key[-1]] = item

        self.parameters['tags'] = tags


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

        if self.parameters['use_opt_wf']:
            use_opt_wf(self.parameters['vmc_in'])



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
