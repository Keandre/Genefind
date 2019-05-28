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
        # Test that arrays of length n <= 2 work
        array = [12]
        array_intact = array.copy()
        random_add_n(array,n=2)
        self.assertTrue(len(array) == len(array_intact) + 2)

        array = [12,16]
        array_intact = array.copy()
        random_add_n(array, n=2)
        self.assertTrue(len(array) == len(array_intact) + 2)
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
        for x in range(100):
            # Test that array mutates
            array = random.sample(range(0, 64), random.randrange(3, 8))
            array_intact = array.copy()
            random_swap(array)
            # Test that the array was altered, and we didn't end up swapping the same element with itself.
            self.assertTrue(array != array_intact)

            # Test that array swaps two cells
            swappedIndexes = []
            for i in range(0, len(array)):
                if array[i] != array_intact[i]:
                    swappedIndexes.insert(0, i)
            self.assertTrue(len(swappedIndexes) == 2)
            self.assertEqual(array[swappedIndexes[0]], array_intact[swappedIndexes[1]])
            self.assertEqual(array[swappedIndexes[1]], array_intact[swappedIndexes[0]])

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
        self.assertGreater(Mod3Run(candidate_a, fitness_a), 0)

        # Test that candidate arrays shorter than fitness array receive score
        # if similar
        candidate_a = [0,1,20]
        self.assertGreater(Mod3Run(candidate_a, fitness_a), 0)

        # Test that candidate array 3/5 of fitness receives score of 60
        candidate_a = [0,1,20]
        fitness_a = [0,1,20,50,30]
        self.assertEqual(Mod3Run(candidate_a, fitness_a), 60)

        # Test that candidate array 5/3 of fitness receives score of 60
        candidate_a = [0,1,20,50,30]
        fitness_a = [0,1,20]
        self.assertEqual(Mod3Run(candidate_a, fitness_a), 60)

        # Test that candidate array 100% of fitness receives score of 100
        candidate_a = [20, 43, 10, 45, 22, 39, 4, 1, 6, 7, 0, 55]
        fitness_a = [0, 55, 20, 43, 10, 45, 22, 39, 4, 1, 6, 7]
        self.assertEqual(Mod3Run(candidate_a, fitness_a), 100)

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


class TestGeneFind(unittest.TestCase):
    """Test functions and classes which run the genefind algorithm."""

    def test_run_genefind(self):
        start_time = time.time()
        print_median_generations_needed(10)
        print("Program took ", str(time.time() - start_time), " time to run")
        # Test that string 'tgattacaa' converges
        self.assertTrue(run_genefind('tgattacaa', 90, 100, 25))
        
class TestParameterFind(unittest.TestCase):
    """Test functions and classes related to the meta-algorithm for finding the best
    parameters to run genefind with."""
    
    def test_fp_cull(self):
        parameter_population = [{"performance":i} for i in range(0,10)]
        # Test that 90th percentile only returns 90th percentile performance
        self.assertEqual(len(fp_cull(parameter_population, .10)),1)
        # Test that algorithm works to 2 digits of precision
        parameter_population = [{"performance":i} for i in range (0,100)]
        self.assertEqual(len(fp_cull(parameter_population,.05)),5)
        
if __name__ == '__main__':
    unittest.main()
