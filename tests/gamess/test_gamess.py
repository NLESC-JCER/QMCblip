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

    @found_gamess
    @found_champ
    def test_userscrError(self):
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        wf = WavefunctionCreator(atoms, str(Path.home().joinpath('software/champ')))
        with self.assertRaises(RuntimeError):
            wf.setup_rhf(userscr=None)
        os.chdir('..')
        with self.assertRaises(FileNotFoundError):
            wf.setup_rhf(userscr="not_exist")
        os.chdir('..')


    @found_champ
    @found_gamess
    def test_gamess(self):
        atoms = Atoms('C2', [(0,0,-0.61385), (0,0,0.61385)])
        wf = WavefunctionCreator(atoms, str(Path.home().joinpath('software/champ')))
        wf.setup_rhf(userscr=str(Path.home().joinpath(Path('software/gamess/restart'))))
        wf.setup_cas(system=dict(mwords=500), drt=dict(nmcc=2, ndoc=2, nval=2), userscr=str(Path.home().joinpath(Path('software/gamess/restart'))))
        wf.setup_ci(system=dict(mwords=500), cidrt=dict(nfzc=2, ndoc=2, nval=2), userscr=str(Path.home().joinpath(Path('software/gamess/restart'))))
        wf.convert_to_champ()
        input = wf.create_champ_input()
        input.optwf.nopt_iter = 10
        input.write('vmc.inp')
        atoms.calc = CHAMP(champ_loc=str(Path.home().joinpath('software/champ'))+"/bin/vmc.mov1", settings = input)
        energy = atoms.get_total_energy()
        print(energy)
        self.assertEqual(int(energy+297.5), 0)

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
