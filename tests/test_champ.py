import os
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
    not Path.home().joinpath(Path('software/champ')).is_dir(), reason="CHAMP not found."
)
class TestChamp(unittest.TestCase):
    def setUp(self):
        self.champ_dir = str(Path.home().joinpath('software/champ'))
        self.settings = Settings.read("tests/test_data/C2_champ/vmc.inp")
        self.atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        os.chdir('tests/test_data')

    def test_setupError(self):
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="")
        with self.assertRaises(CalculatorSetupError):
            t = CHAMP(champ_loc=".vmc.mov1")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", vmc_in="C2_champ/vmc.inp")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", settings=self.settings)
            remove('vmc.inp')

    def test_writeInput(self):
        calc = CHAMP(champ_loc="./vmc.mov1", settings=self.settings)
        calc.write_input(atoms=self.atoms)
        f = open('molecule.xyz')
        self.assertEqual(int(f.readline()), len(self.atoms))
        f.readline()
        for i in range(len(self.atoms)):
            a = np.array(f.readline().split()[1:], dtype=float) * units.Bohr/units.Ang
            b = self.atoms[i].position
            for j in range(3):
                self.assertAlmostEqual(a[j], b[j])
        remove('vmc.inp')
        remove('vmc_temp.inp')
        remove('molecule.xyz')
        f.close()

    @found_champ
    def test_C2(self):
        calc = CHAMP(champ_loc=self.champ_dir+"/bin/vmc.mov1", settings=self.settings)
        self.atoms.calc = calc

        self.assertAlmostEqual(self.atoms.get_total_energy(), -292.4598135740918)
        cleanup()

    def tearDown(self):
        os.chdir("../..")



if __name__ == '__main__':
    unittest.main()
