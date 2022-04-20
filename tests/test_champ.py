import os
import shutil
import unittest
from os import remove
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.calculator import CalculatorSetupError
from ase.md.verlet import VelocityVerlet
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings, cleanup

found_champ = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/champ')).is_dir(), reason="CHAMP not found."
)
class TestChamp(unittest.TestCase):
    def setUp(self):
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')

    def test_setupError(self):
        settings = Settings.read("../C2_champ/vmc.inp")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="")
        with self.assertRaises(CalculatorSetupError):
            t = CHAMP(champ_loc="../vmc.mov1")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", vmc_in="../C2_champ/vmc.inp")
        with self.assertRaises(CalculatorSetupError):
            CHAMP(champ_loc="", settings=settings)

    def test_writeInput(self):
        settings = Settings.read("../C2_champ/vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        calc = CHAMP(champ_loc="../vmc.mov1", settings=settings)
        calc.write_input(atoms=atoms)
        f = open('molecule.xyz')
        self.assertEqual(int(f.readline()), len(atoms))
        f.readline()
        for i in range(len(atoms)):
            a = np.array(f.readline().split()[1:], dtype=float) * units.Bohr/units.Ang
            b = atoms[i].position
            for j in range(3):
                self.assertAlmostEqual(a[j], b[j])
        f.close()

    @found_champ
    def test_C2(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")
        atoms.calc = calc

        self.assertAlmostEqual(atoms.get_total_energy(), -292.4598135740918)

    @found_champ
    def test_C2_MD(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")
        dyn = VelocityVerlet(atoms, units.fs)
        dyn.run(3)
        self.assertListEqual(atoms.get_positions().tolist, [[0,0,-0.61],[0,0,0.61]])
        self.assertAlmostEqual(atoms.get_total_energy(), -292.4598135740918)

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
