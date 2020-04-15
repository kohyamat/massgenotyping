import unittest

from massgenotyping.find_ssrs import rotate


class TestRotate(unittest.TestCase):
    def test_rotate(self):
        ans = ["test", "estt", "stte", "ttes"]
        for x, y in zip(rotate("test"), ans):
            self.assertEqual(
                x, y,
            )


if __name__ == "__main__":
    unittest.main()
