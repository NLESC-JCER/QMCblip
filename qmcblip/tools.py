import numpy as np

from flare import otf_parser

import matplotlib.pyplot as plt

from ase.db import connect
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.calculators.calculator import Calculator

import periodictable as pt

from os.path import exists
from os import remove


class Analyze():
    """Tool for analyzing OTF data.

    Args:
        file (Path): otf file to analyze.

    Attributes:
        results (dict): dictionary containing simulation results.
    """

    def __init__(self, file):
        self.otf = otf_parser.OtfAnalysis(file, calculate_energy=True)
        self.results = None

    def to_xyz(self, filename=""):
        """Create an .xyz format file from the OTF trajectory.

        Args:
            filename (str): filename for the .xyz file (include extension).
        """
        if (filename == ""):
            self.otf.to_xyz('traj.xyz')
        else:
            self.otf.to_xyz(filename)

    def get_data(self):
        """Retrieve the data from the OTF file.

        Note:
            The data will be stores in the results attributes.
            Allowed keys are: times, potential energy, kinetic energy, total energy, temperature.
        """
        dt = self.otf.header['dt']

        frames = np.arange(len(self.otf.times)+1)

        self.results = dict()

        self.dft = self.otf.dft_frames
        self.nondft = [i for i in frames if i not in self.dft]

        totE = np.zeros(len(frames))
        kinE = np.zeros(len(frames))
        potE = np.zeros(len(frames))
        temp = np.zeros(len(frames))

        dftframe = 0
        for i in range(len(frames)):
            if i in self.dft:
                totE[i] = self.otf.gp_thermostat['total energy'][dftframe]
                kinE[i] = self.otf.gp_thermostat['kinetic energy'][dftframe]
                potE[i] = self.otf.gp_thermostat['potential energy'][dftframe]
                temp[i] = self.otf.gp_thermostat['temperature'][dftframe]
                dftframe += 1
            else:
                totE[i] = self.otf.thermostat['total energy'][i-1]
                kinE[i] = self.otf.thermostat['kinetic energy'][i-1]
                potE[i] = self.otf.thermostat['potential energy'][i-1]
                temp[i] = self.otf.thermostat['temperature'][i-1]


        self.results['times'] = frames[0:-2] * dt
        self.results['potential energy'] = potE[1:-1]
        self.results['kinetic energy'] = kinE[2:]
        self.results['total energy'] = kinE[2:] + potE[1:-1]
        self.results['temperature'] = temp[2:]

        if (max(self.dft) > max(self.nondft)):
            self.dft = self.dft[:-1]
        else:
            self.nondft = self.nondft[:-1]

        return self.results
    
    def plot_energy(self, filename=""):
        """Plot the energy

        Args:
            filename (str): file to save the energy plot to (leave empty if not wanted).
        """
        if self.results is None:
            self.get_data()

        times = self.results['times']
        totE = self.results['total energy']
        potE = self.results['potential energy']
        kinE = self.results['kinetic energy']
        plt.figure()
        plt.grid()
        plt.title("Energy")
        plt.plot(times, totE - totE[0], label="Total Energy")
        plt.plot(times, potE - potE[0], label="Potential Energy")
        plt.plot(times, kinE - kinE[0], label="Kinetic Energy")
        plt.scatter([times[index-1] for index in self.dft[1:]], [totE[index-1] - totE[0] for index in self.dft[1:]], marker='x', label="QMC Calculations", color='black')
        plt.legend()
        plt.xlabel("Time (ps)")
        plt.ylabel("Energy (eV)")
        if (filename != ""):
            plt.savefig(filename)

    def plot_temperature(self, filename=""):
        """Plot the temperature

        Args:
            filename (str): file to save the temperature plot to (leave empty if not wanted).
        """
        if self.results is None:
            self.get_data()

            plt.figure()

        times = self.results['times']
        temp = self.results['temperature']
        plt.figure()
        plt.grid()
        plt.title("Temperature")
        plt.plot(times, temp)
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")
        if (filename != ""):
            plt.savefig(filename)

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
    db = connect(db_file)
    traj = Trajectory(traj_file)

    # Write to database
    for atoms in traj:
        db.write(atoms, relaxed=False)

def otf_to_db(otf_file, db_file, append=False):
    # File checks
    if not exists(otf_file):
        raise FileNotFoundError(otf_file + " is not in this directory!")
    if exists(db_file) and append:
        print(db_file + " already exists. Appending to it...")
    elif exists(db_file) and not append:
        print(db_file + " already exists. Overwriting...")
        remove(db_file)

    # Open files
    db = connect(db_file)
    otf = otf_parser.OtfAnalysis(otf_file, calculate_energy=True)

    types = otf.header['species']
    atomicnumbers = np.zeros(len(types))
    atomicmasses = np.zeros(len(types))
    for i in range(len(types)):
        atomicnumbers[i] = getattr(pt, types[i]).number
        atomicmasses[i] = getattr(pt, types[i]).mass

    for i in range(len(otf.dft_frames)):
        forces = otf.gp_force_list[i]
        positions = otf.gp_position_list[i]
        energy = otf.gp_thermostat['potential energy'][i]

        atoms = Atoms(''.join(types), positions=positions)
        atoms.calc = Calculator()
        atoms.calc.implemented_properties = ['forces', 'energy']
        atoms.calc.results['forces'] = forces
        atoms.calc.results['energy'] = energy
        atoms.numbers = atomicnumbers
        atoms.masses = atomicmasses
        
        db.write(atoms)