import os
import shutil
import unittest
from pathlib import Path

import pytest
from ase import Atoms
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings
from qmcblip.gamess.utils import WavefunctionCreator

found_champ = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/champ')).is_dir(), reason="CHAMP not found."
)
found_games = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/gamess')).is_dir(), reason="GAMESS not found."
)

class TestGamess(unittest.TestCase):

    def setUp(self):
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')

    @found_champ
    @found_games
    def test_gamess(self):
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        wf = WavefunctionCreator(atoms, str(Path.home().joinpath('software/champ')))
        wf.setup_rhf()
        wf.setup_cas(system=dict(mwords=500), drt=dict(nmcc=2, ndoc=2, nval=2))
        wf.setup_ci(system=dict(mwords=500), cidrt=dict(nfzc=2, ndoc=2, nval=2))
        wf.convert_to_champ()
        input = wf.create_champ_input()
        input.optwf.nopt_iter = 10
        input.write('vmc.inp')
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1", settings = input)

        self.assertAlmostEqual(atoms.get_total_energy(), -297.5, places=0)

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
