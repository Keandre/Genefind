import unittest
from genefind import *
import pdb

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

    def test_random_swap(self):
        # Test that array mutates
        array = [10,12,45,62,4]
        array_intact = array.copy()
        random_swap(array)
        # TODO: Interesting bug here where occasionally same array returned
        self.assertTrue(array != array_intact)

        # Test that array swaps two cells

    def test_run_2_mod(self):
        # Make sure that the population is not literally all the same array object
        # -__-
        output = Mod2Run([[12,10,20]],50)
        self.assertNotEqual(id(output[0]),id(output[1]))
        
    def test_run_3_mod(self):
        # Test that candidate arrays longer than fitness array receive score if
        # similar
        candidate_a = [0,1,20,45,52,10,8,4,3,1]
        fitness_a = [0,1,20,50,30,40,10,60]
        self.assertTrue(Mod3Run(candidate_a, fitness_a) > 0)

        # Test that candidate arrays shorter than fitness array receive score
        # if similar
        candidate_a = [0,1,20]
        self.assertTrue(Mod3Run(candidate_a, fitness_a) > 0)

        # Test that candidate array 3/5 of fitness receives score of 60
        candidate_a = [0,1,20]
        fitness_a = [0,1,20,50,30]
        self.assertTrue(Mod3Run(candidate_a, fitness_a) == 60)

        # Test that candidate array 5/3 of fitness receives score of 60
        candidate_a = [0,1,20,50,30]
        fitness_a = [0,1,20]
        self.assertTrue(Mod3Run(candidate_a, fitness_a) == 60)

        # Test that candidate array 100% of fitness receives score of 100
        candidate_a = [0, 55, 20, 43, 10, 45, 22, 39, 4, 1, 6, 7]
        fitness_a = [0, 55, 20, 43, 10, 45, 22, 39, 4, 1, 6, 7]
        self.assertTrue(Mod3Run(candidate_a, fitness_a) == 100)

    def test_genetic2b64(self):
        codons = {"AAA":0, "AAT":1, "AAG":2,
                  "AAC":3, "ATA":4, "ATT":5,
                  "ATG":6, "ATC":7, "AGA":8,
                  "AGT":9, "AGG":10, "AGC":11,
                  "ACA":12, "ACT":13, "ACG":14,
                  "ACC":15, "TAA":16, "TAT":17,
                  "TAG":18, "TAC":19, "TTA":20,
                  "TTT":21, "TTG":22, "TTC":23,
                  "TGA":24, "TGT":25, "TGG":26,
                  "TGC":27, "TCA":28, "TCT":29,
                  "TCG":30, "TCC":31, "GAA":32,
                  "GAT":33, "GAG":34, "GAC":35,
                  "GTA":36, "GTT":37, "GTG":38,
                  "GTC":39, "GGA":40, "GGT":41,
                  "GGG":42, "GGC":43, "GCA":44,
                  "GCT":45, "GCG":46, "GCC":47,
                  "CAA":48, "CAT":49, "CAG":50,
                  "CAC":51, "CTA":52, "CTT":53,
                  "CTG":54, "CTC":55, "CGA":56,
                  "CGT":57, "CGG":58, "CGC":59,
                  "CCA":60, "CCT":61, "CCG":62,
                  "CCC":63}
        # Test that string is serialized to proper length
        genstring = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        array_list = genetic2b64(genstring, codons)
        self.assertEqual(len(array_list), len(genstring) / 3)
        
if __name__ == '__main__':
    unittest.main()
