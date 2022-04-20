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
    not Path('~/software/champ').is_dir()
)
class TestChamp(unittest.TestCase):
    def setUp(self):
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
        calc = CHAMP(champ_loc=".vmc.mov1", settings=self.settings)
        calc.write_input(atoms=self.atoms)
        f = open('molecule.xyz')
        self.assertEqual(int(f.readline()), len(self.atoms))
        f.readline()
        for i in range(len(self.atoms)):
            a = np.array(f.readline().split()[1:], dtype=float) * units.Bohr/units.Ang
            b = self.atoms[i].position
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

        self.assertAlmostEqual(atoms.get_total_energy(), -293.0584640666279)
        cleanup()




if __name__ == '__main__':
    unittest.main()
