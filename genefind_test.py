import unittest
from genefind import random_delete_n
from genefind import random_add_n

#TODO: Add class testing initial and other misc functions

class TestMutations(unittest.TestCase):

    def test_random_delete_n(self):
        # Test that arrays longer than n+1 are different
        array = [10,12,45,62,4]
        array_intact = array.copy()
        n = 2
        self.assertTrue(array_intact != random_delete_n(array, n=n))
        # Test that n cells were deleted
        array = array_intact.copy()
        self.assertEqual(len(array_intact) - n, len(random_delete_n(array,n=n)))
        # Test consecutive deletion
        array = [10,11,12]
        random_delete_n(array, n=n, consecutive=True)
        #TODO: Currently fails because index meanings change as array size changes
        self.assertTrue(array == [10] or array == [12])

    def test_random_add_n(self):
        # Test that non-consecutive adds work
        array = [10,12,45,62,4]
        array_intact = array.copy()
        random_add_n(array, n=3)
        self.assertTrue(len(array) == len(array_intact) + 3)
        # Test that all values are less than 64
        [self.assertTrue(cell < 64) for cell in array]
        # Test that consecutive adds work
        array = [10]
        random_add_n(array, n=3, consecutive=True)
        self.assertTrue(array[0] == 10 or array[-1] == 10)


if __name__ == '__main__':
    unittest.main()
