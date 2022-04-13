import unittest
from qmcblip.champio import Settings
from pydantic import ValidationError

class TestChampio(unittest.TestCase):
    def test_SetupError(self):
        with self.assertRaises(ValidationError):
            Settings()


if __name__ == '__main__':
    unittest.main()