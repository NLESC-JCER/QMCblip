import os
import shutil
import unittest

import pytest
from qmcblip.sgdml.utils import db_to_sgdml


class TestSgdml(unittest.TestCase):

    def setUp(self):
        os.mkdir("tests/test_data/temp")
        os.chdir('tests/test_data/temp')

    def test_db_to_sgdml(self):
        shutil.copyfile("../test.db", "test.db")
        try:
            db_to_sgdml('test.db', 'test.npz')
        except:
            self.fail("Cannot convert DB to dataset.")

    def tearDown(self):
        os.chdir("../../..")
        shutil.rmtree('tests/test_data/temp')


if __name__ == '__main__':
    unittest.main()
