import unittest
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings
from ase.calculators.calculator import CalculatorSetupError
from os import remove
from ase import Atoms, units
import numpy as np


class TestChamp(unittest.TestCase):
    def setUp(self):
        self.settings = Settings.read("examples/C2-quicksim/vmc.inp")
        self.atoms = [Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)]),
        Atoms('C4SH4', [(0,0.71495093597,1.31902341514), (0,-0.71495093597,1.31902341514), (0,-1.24055334534,0.05119099544), (0,1.24055334534,0.05119099544), (0,0,-1.15285178278), (0,1.32194923477,2.21441153704), (0,-1.32194923477,2.21441153704), (0,2.27909548764,-0.24288695123), (0,-2.27909548764,-0.24288695123)])]

    def test_setupError(self):
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="")
        with self.assertRaises(CalculatorSetupError):
            t = CHAMP(champ_loc="tests/vmc.mov1")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", vmc_in="examples/C2-quicksim/vmc.inp")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", settings=self.settings)
            remove('vmc.inp')

    def test_writeInput(self):
        calc = CHAMP(champ_loc="tests/vmc.mov1", vmc_in="examples/C2-quicksim/vmc.inp")
        for atom in self.atoms:
            calc.write_input(atoms=atom)
            f = open('molecule.xyz')
            self.assertEqual(int(f.readline()), len(atom))
            f.readline()
            for i in range(len(atom)):
                a = np.array(f.readline().split()[1:], dtype=float) * units.Bohr/units.Ang
                b = atom[i].position
                for j in range(3):
                    self.assertAlmostEqual(a[j], b[j])
            remove('vmc_temp.inp')
            remove('molecule.xyz')
            f.close()

if __name__ == '__main__':
    unittest.main()