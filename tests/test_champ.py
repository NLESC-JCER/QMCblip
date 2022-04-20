import unittest
from os import remove
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.calculator import CalculatorSetupError
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings, cleanup

found_champ = pytest.mark.skipif(
    not Path('~/software/champ').is_dir()
)
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
    
    @found_champ
    def test_C2(self):
        calc = CHAMP(champ_loc="~/software/champ/bin/vmc.mov1", settings=self.settings)
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        atoms.calc = calc
        atoms.calc.parameters.settings.optwf.vmc_nstep = 2
        atoms.calc.parameters.settings.optwf.vmc_nblk = 2

        self.assertAlmostEqual(atoms.get_total_energy(), -300)
        cleanup()




if __name__ == '__main__':
    unittest.main()
