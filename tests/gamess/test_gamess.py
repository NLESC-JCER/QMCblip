import os
import shutil
import unittest
from pathlib import Path

import pytest
from ase import Atoms
from qmcblip.champ import CHAMP
from qmcblip.gamess.utils import Presets

found_champ = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/champ')).is_dir(), reason="CHAMP not found."
)
found_gamess = pytest.mark.skipif(
    not Path.home().joinpath(Path('software/gamess')).is_dir(), reason="GAMESS not found."
)

class TestGamess(unittest.TestCase):

    def setUp(self):
        if Path.home().joinpath(Path('software/gamess')).is_dir():
            shutil.rmtree(Path.home().joinpath(Path('software/gamess/restart')))
            os.mkdir(Path.home().joinpath(Path('software/gamess/restart')))
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')

    @found_champ
    @found_gamess
    def test_gamess(self):
        atoms = Presets.C2().atoms
        input = Presets.C2(userscr=str(Path.home().joinpath(Path('software/gamess/restart')))).build(str(Path.home().joinpath('software/champ')))
        input.optwf.nopt_iter = 10
        input.write('vmc.inp')
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1", settings = input)
        energy = atoms.get_total_energy()
        print(energy)
        self.assertEqual(int((energy+297.5)/5), 0)

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
