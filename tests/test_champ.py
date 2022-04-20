import os
import shutil
import unittest
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.calculator import CalculatorSetupError
from ase.md.verlet import VelocityVerlet
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings

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
    def test_C2_2n(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1", ncore=2)
        atoms.calc = calc
        print(atoms.get_total_energy())
        self.assertAlmostEqual(atoms.get_total_energy(), -292.5087940689357)

    @found_champ
    def test_C2_MD(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")
        dyn = VelocityVerlet(atoms, units.fs)
        dyn.run(3)
        res = [[0.004348154138868875, 0.0021912427163622832, -0.5825713579980879], 
               [-0.0025218801616919595, 0.001956434734404338, 0.5759901971368933]]
        for i in range(2):
            for j in range(3):
                self.assertAlmostEqual(atoms.get_positions()[i][j], res[i][j])
        self.assertAlmostEqual(atoms.get_total_energy(),  -293.142893130546)

    @found_champ
    def test_C2_MD_opt_wf(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1", use_opt_wf=True)
        dyn = VelocityVerlet(atoms, units.fs)
        dyn.run(3)
        res = [[ 6.57243666e-04, -1.50591101e-03, -5.79840171e-01],
               [ 3.53073963e-04, -2.79268146e-03, 5.80160063e-01]]
        print(atoms.get_positions(), atoms.get_total_energy())
        for i in range(2):
            for j in range(3):
                self.assertAlmostEqual(atoms.get_positions()[i][j], res[i][j])
        self.assertAlmostEqual(atoms.get_total_energy(),  -295.01571168245135)

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
