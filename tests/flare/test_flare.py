import os
import shutil
import unittest
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units
from ase.atoms import Cell
from qmcblip.champ import CHAMP
from qmcblip.flare.quicksim import OTFSettings, quicksim
from qmcblip.flare.utils import Analyze

found_champ = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/champ')).is_dir(), reason="CHAMP not found."
)

class TestChamp(unittest.TestCase):

    def setUp(self):
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')
        np.random.seed(123)

    @found_champ
    def test_quicksim_C2(self):
        shutil.copytree("../C2_champ/pool", "pool")
        shutil.copyfile("../C2_champ/vmc.inp", "vmc.inp")

        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
        atoms.pbc=[True, True, True]

        settings = OTFSettings(theory=OTFSettings.FLARE())
        settings.std_tolerance_factor = 3

        calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")

        quicksim(atoms, 0.5, 5, calc, settings)

        res = Analyze('OTF.out')
        res.to_xyz()
        res.get_data()
        self.assertAlmostEqual(res.results['total energy'][-1], -292.466862)

    @found_champ
    def test_quicksim_Thio(self):
        shutil.copytree("../Thio_champ/pool", "pool")
        shutil.copyfile("../Thio_champ/vmc.inp", "vmc.inp")

        atoms = Atoms('C4SH4', [(0,0.71495093597,1.31902341514), (0,-0.71495093597,1.31902341514),
            (0,-1.24055334534,0.05119099544), (0,1.24055334534,0.05119099544),
            (0.6456,0,-1.15285178278), (0,1.32194923477,2.21441153704), (0,-1.32194923477,2.21441153704),
            (0,2.27909548764,-0.24288695123), (0,-2.27909548764,-0.24288695123)])
        atoms.cell  = Cell.fromcellpar([50, 50, 50, 90, 90, 90])
        atoms.pbc=[True, True, True]

        settings = OTFSettings(theory=OTFSettings.FLAREPP())
        settings.std_tolerance_factor = 3

        calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")

        quicksim(atoms, 0.5, 5, calc, settings)

        res = Analyze('OTF.out')
        res.to_xyz()
        res.get_data()
        self.assertAlmostEqual(res.results['total energy'][-1], -956.934665)



    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
