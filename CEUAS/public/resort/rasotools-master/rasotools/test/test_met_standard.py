import unittest

import rasotools


class Standard(unittest.TestCase):
    def test_default(self):
        self.assertFalse(False)


if __name__ == '__main__':
    unittest.main()
