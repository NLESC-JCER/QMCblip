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
        settings.std_tolerance_factor = 0.5

        calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1")

        quicksim(atoms, 0.5, 5, calc, settings)

        res = Analyze('OTF.out')
        res.to_xyz()
        res.get_data()
        print(res.results['total energy'][-1])



    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
