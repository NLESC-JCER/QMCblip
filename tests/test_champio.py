import os
import shutil
import unittest

from pydantic import ValidationError
from qmcblip.champio import Settings


class TestChampio(unittest.TestCase):
    def setUp(self):
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')

    def test_SetupError(self):
        with self.assertRaises(ValidationError):
            Settings()

    def test_write(self):
        settings = Settings.read('../C2_champ/vmc.inp')
        settings.extra_obj = 'test'
        try:
            settings.write('temp.inp')
        except:
            self.fail("Cannot write vmc input with custom tags.")

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
