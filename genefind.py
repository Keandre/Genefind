#mutations:
#dshift by [val={1,3}]
#drandom delete 1
#drandom delete 3 consecutive
#drandom delete 3 random (essentially RD1 three times)
#drandom add 1
#drandom add 3 consecutive
#drandom add 3 random
#drandom switch

from argparse import ArgumentParser
import multiprocessing
import random
import numpy as np
import statistics
from functools import partial
import json
import pprint
import cProfile
import time
import pdb

#interesting functions:
#random.randrange()
#random.choice()
#random.choices()
#random.sample()

#algo
#input = a list of strings of ACTG seq
#output = a list of a list of 100 strings of mutated ACTG seq
#If the total number of mutations applied to the total number of input strings
#is less than 100, begin using output strings as input strings, in the order
#they were written.

actg = "ACTG"
ACID_SPACE = 64

def random_delete_n(array, n=1, consecutive=False):
    """Randomly delete n cells from a array."""
    if len(array)<n+1: return array
    if not consecutive:
        for i in range(n):
            cell_i = random.randrange(len(array))
            del(array[cell_i])
        return array
    else:
        start_cell_i = random.randrange(len(array)-n)
        cells_deleted = 0 
        for i in range(start_cell_i,start_cell_i + n):
            del(array[i - cells_deleted])
            cells_deleted += 1
        return array

def random_add_n(array, n=1, consecutive=False):
    """Randomly add n cells to a array."""
    if not consecutive:
        for i in range(n):
            insert_i = random.randrange(len(array))
            array.insert(insert_i,
                         random.randrange(ACID_SPACE))
        return array
    else:
        insert_i = random.randrange(len(array))
        for i in range(n):
            array.insert(insert_i,
                         random.randrange(ACID_SPACE))
        return array

def random_swap(array):
    """Randomly swap the contents of two cells."""
    if len(array)<2: return array
    idx = random.sample(range(len(array)), 2)
    array[idx[0]], array[idx[1]] = array[idx[1]], array[idx[0]]
    return array

def random_shift(stringIn, n):
    """Deterministically shift a string starting from position n.

    It's important to note that how this function is used in the rest of the program
    relies on you specifying constant values for the random shifts, so in the 
    set of possible random operations positions are defined such that it can be 
    1 or 3, positive or negative. This allows it to implement shifting both sides
    of the string in both lengths with one function."""
    if len(stringIn)<n+1: return stringIn
    stringOut = stringIn[n:]+stringIn[:n]
    return stringOut

# Gotten from https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/ 
def lcs(X , Y): 
    """Dynamic Programming implementation of Longest Common Substring problem."""
    # find the length of the strings 
    m = len(X)
    n = len(Y)
  
    # declaring the array for storing the dp values 
    L = [[None]*(n+1) for i in range(m+1)]
  
    """Following steps build L[m+1][n+1] in bottom up fashion 
    Note: L[i][j] contains length of LCS of X[0..i-1] 
    and Y[0..j-1]"""
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1]+1
            else: 
                L[i][j] = max(L[i-1][j], L[i][j-1])
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return L[m][n] 

# Here we define a set of partial functions that are used as operations to mutate
# strings. 
operations = [partial(random_shift,n=3), partial(random_shift,n=1),
              partial(random_shift,n=-1),
              partial(random_shift,n=-3), random_swap,
              random_add_n, partial(random_add_n,n=3),
              partial(random_add_n,n=3,consecutive=True),
              random_delete_n,partial(random_delete_n,n=3),
              partial(random_delete_n, n=3,consecutive=True),]

def Mod2Run(listIn, population):
    """Generate the population. Here we randomly pull from a set of predefined
    operations to mutate our existing strings, then execute the operation."""
    countProc = 0
    listOut = []
    listOut += listIn
    while len(listOut)<=population:
        array = listOut[countProc].copy()
        if len(array) < 2:
            # If we've evolved to extinction, add new stuff to select on
            array[0] = random.randrange(ACID_SPACE)
            array = random_add_n(array, n=2)
            listOut.append(array)
            countProc += 1
        else:
            operation = random.choice(operations)
            listOut.append(operation(array))
            countProc+=1
    return listOut

def Mod3Run(candidate_a, fitness_a):
    """Test fitness by comparing the similarity of two arrays.

    The candidate array is compared to the fitness array, and a fitness score is 
    returned by the procedure."""
    candidate_length = len(candidate_a)
    fitness_length = len(fitness_a)
    comparison_length = max(candidate_length, fitness_length)
    match = lcs(candidate_a, fitness_a)
    return (match/comparison_length)*100

def Mod4Run(listIn, fitness, fitness_percentile):##pop cull
    listP = []
    for i in listIn:
        listP.append([Mod3Run(i,fitness),i])
    listP.sort()
    fitness_cutoff = -round(len(listP) * fitness_percentile)
    if not fitness_cutoff: # Handle case where we round to zero
        fitness_cutoff = 1
    return [candidate for candidate in reversed(listP[fitness_cutoff:])]

def Mod1Run(listIn, trackIn={}):##iterator
    trackOut = trackIn
    bestE = listIn[0]
    listOut = []
    if len(listIn) == 1:
        return listIn, trackOut
    for i in listIn:
        listOut.append(i[1])
    
    trackOut["lastFit"] = trackOut["currentFit"]
    trackOut["currentFit"] = bestE[0]
    trackOut["generationCount"]+=1
    bestFit = trackOut["bestFit"]
    if bestFit<bestE[0]:
        trackOut["bestFit"]=bestE[0]
    if bestE[0] == 100:
        print("Result Found.","\n"\
              ,"Generation count =", trackOut["generationCount"], '\n'\
              ,"String :\n", bestE[1], GetGeneString(bestE[1]))
        return listOut, trackOut
    return listOut, trackOut

def GetGeneString(codon_list):
    gene_string = ""
    codon_keys = list(codons.keys())
    for c in codon_list:
        gene_string += codon_keys[c]
    return gene_string

def initial(initial_array_size):
    array = []
    for i in range(initial_array_size):
        array.append(
            random.randrange(ACID_SPACE)
            )
    return array

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

def genetic2b64(genstring, codons):
    return [codons[genstring[int(3 * (i/3 - 1)):i]] for i in range(3,len(genstring) + 3,3)]

##base program
print("Modules loaded.")
print(\
    """
********************
"GeneFind: A Designspace Searcher"
********************
""")


def run_genefind(target_string, fitness_percentile, population, initial_string_size):
    """Run one instance of the genefind algorithm and return the number of 
    generations to find the target string."""
    lis = [initial(initial_string_size)]
    fitness = genetic2b64(target_string.upper(), codons)
    fitness_cutoff = (100 - fitness_percentile) / 100
    track = {"generationCount":0, 'currentFit':0.0, \
             'lastFit':0.0, 'bestFit':0.0}
    while 1:
        lis, track = Mod1Run(lis, track)
        lis = Mod2Run(lis, population)
        lis = Mod4Run(lis,fitness, fitness_cutoff)
        if track["bestFit"] == 100:
            break
        if track["generationCount"]%1000 == 0 or track["bestFit"]>track["lastFit"]:
            #TODO: This should really be put under a -v option
            print("Generation ",track["generationCount"],"Fitness ", track["bestFit"])
    return track["generationCount"] * population

class ThreadableGenefind:
    """An OOP version of the genefind mainloop which exposes its important tracking 
    values as attributes."""
    def run(self, target_array, fitness_percentile,
            population, initial_string_size, shared_state=None):
        if shared_state:
            shared_threads = shared_state[0]
            thread_id = shared_state[1]
            shared_threads[thread_id] = 0
            out_queue = shared_state[2]
            self.ready = shared_state[3]
            self.ready.acquire()
        gene_population = [initial(initial_string_size)]
        fitness = target_array
        fitness_cutoff = (100 - fitness_percentile) / 100
        self.population = population
        self.generation_count = 0
        self.current_fit = 0.0
        self.last_fit = 0.0
        self.best_fit = 0.0
        self._stopped = False
        while 1:
            # Update external performance metric
            if shared_threads[thread_id] is False:
                self._stopped = True
            else:
                shared_threads[thread_id] = self.performance()
            if self._stopped:
                break
            if self.best_fit == 100:
                break
            self.generation_count += 1
            gene_population = self.record(gene_population)
            gene_population = Mod2Run(gene_population, population)
            gene_population = Mod4Run(gene_population, fitness, fitness_cutoff)

            #if self.generation_count % 1000 == 0 or self.best_fit > self.last_fit:
                # This should really be put under a -v option
                #print("Generation ", self.generation_count, "Fitness ", self.best_fit)
        self._stopped = True
        try:
            out_queue.put(self.performance())
            self.ready.release()
        except NameError:
            return self.performance()

    def record(self, gene_population):
        """Do record keeping on the main loop of the algorithm. Update attributes
        and check to see if we've found the result we're looking for."""
        bestE = gene_population[0]
        if len(gene_population) == 1:
            return gene_population
        gene_population_out = [gene[1] for gene in gene_population]
        self.last_fit = self.current_fit
        self.current_fit = bestE[0]
        bestFit = self.best_fit
        if bestFit<bestE[0]:
            self.best_fit = bestE[0]
            if bestE[0] == 100:
                print("Result Found.","\n"\
                      ,"Generation count =", self.generation_count, '\n'\
                      ,"String :\n", bestE[1])
        return gene_population_out

    def performance(self):
        return self.generation_count * self.population

    def stopped(self):
        return self._stopped
    
    def stop(self):
        self._stopped = True

def fp_cull(parameter_population, fitness_percentile):
    """Cull the parameter population to find its top members above the 
    fitness percentile."""
    parameter_population.sort(key=(lambda candidate: candidate["performance"]))
    fitness_cutoff = -round(len(parameter_population) * fitness_percentile)
    if not fitness_cutoff: # Handle case where we round to zero
        fitness_cutoff = 1
    return [candidate for candidate in reversed(parameter_population[fitness_cutoff:])]
    
def fp_sample(target_string, parameter_population, sample_runs,
              performance_limit, num_threads):
    """Sample the parameter population in the find-parameters routine by running
    each candidate parameter set multiple times and taking the median performance
    as its score."""
    thread_count = sample_runs
    for candidate in parameter_population:
        sample_runs = thread_count
        trials = []
        trials_queue = multiprocessing.Queue()
        threads_host = {}
        manager = multiprocessing.Manager()
        threads_shared = manager.dict()
        race_controller = manager.BoundedSemaphore(value=num_threads)
        for run in range(sample_runs):
            sample = ThreadableGenefind()
            run = multiprocessing.Process(target=sample.run,
                                          args=(target_string,
                                                candidate["fitness_percentile"],
                                                candidate["population"],
                                                candidate["string_size"],
                                                (threads_shared,
                                                 sample_runs,
                                                 trials_queue,
                                                 race_controller)))
            run.start()
            threads_host[sample_runs] = run
        while sample_runs:
            time.sleep(0.1)
            try:
                trials.append(trials_queue.get(block=False))
                sample_runs -= 1
            except:
                None
            for thread in threads_shared.keys():
                if threads_shared[thread] > performance_limit:
                    #TODO: Make this not an ugly hack
                    threads_shared[thread] = False # Send message to stop the thread
        candidate["performance"] = statistics.median(trials)
        for thread in threads_host.values():
            thread.terminate()
    return parameter_population
        

def fp_mutate(parameter_population, population):
    """Increase the parameter population to its full size, then mutate it."""
    original_population = parameter_population.copy()
    if len(parameter_population) < population:
        for index in range(population - len(parameter_population)):
            parameter_population.append(random.choice(original_population).copy())
    for candidate in parameter_population:
        # Mutate fitness percentile
        sign = random.randrange(2)
        mutation = random.randint(1,9) * 0.1
        if sign:
            if candidate["fitness_percentile"] + mutation > 99.9:
                candidate["fitness_percentile"] -= mutation
            else:
                candidate["fitness_percentile"] += mutation
        else:
            if candidate["fitness_percentile"] - mutation < 0.1:
                candidate["fitness_percentile"] += mutation
            else:
                candidate["fitness_percentile"] -= mutation

        # Mutate population
        sign = random.randrange(2)
        mutation = random.randint(1,9)

        if sign:
            if candidate["population"] + mutation > 1000:
                candidate["population"] -= mutation
            else:
                candidate["population"] += mutation
        else:
            if candidate["population"] - mutation < 1:
                candidate["population"] += mutation
            else:
                candidate["population"] -= mutation

        # Mutate string_size
        sign = random.randrange(2)
        mutation = random.randint(1,9)

        if sign:
            if candidate["string_size"] + mutation > 1000:
                candidate["string_size"] -= mutation
            else:
                candidate["string_size"] += mutation
        else:
            if candidate["string_size"] - mutation < 1:
                candidate["string_size"] += mutation
            else:
                candidate["string_size"] -= mutation
        
    return parameter_population
            
            

def fp_record(tracking_dict, parameter_population):
    tracking_dict["last_best"] = tracking_dict["current_best"]
    tracking_dict["current_best"] = parameter_population[-1]
    if (tracking_dict["absolute_best"]["performance"] <
        tracking_dict["current_best"]["performance"]):
        tracking_dict["absolute_best"] = parameter_population[-1]
    tracking_dict["history"].append(parameter_population)
    return tracking_dict

def find_parameters(target_string, generations=10, sample_runs=5,
                    population=5, fitness_percentile=90, num_threads=1,
                    outpath=None):
    """Find the parameters that minimize the number of individuals for genefind.

    This is done using a genetics algorithm that starts with an initial parameter
    set, then mutates and finds the fitness of each member, finally culling the 
    members which do not meet a certain level of relative performance."""
    fitness_cutoff = (100 - fitness_percentile) / 100
    fitness = genetic2b64(target_string.upper(), codons)
    initial_parameters = {"fitness_percentile":90,
                          "population":100,
                          "string_size":50,
                          "performance":518600}
    parameter_population = [initial_parameters]
    tracking = {"generation_count":0,
                "absolute_best":initial_parameters,
                "current_best":initial_parameters,
                "last_best":initial_parameters,
                "history":[]}
    for generation in range(generations):
        tracking["generation_count"] += 1
        tracking = fp_record(tracking, parameter_population)
        print("Generation {}: Current best fit is ({},{},{}) with performance {}.".format(
            tracking["generation_count"],
            tracking["current_best"]["fitness_percentile"],
            tracking["current_best"]["population"],
            tracking["current_best"]["string_size"],
            tracking["current_best"]["performance"]))
        parameter_population = fp_mutate(parameter_population, population)
        performance_limit = tracking["current_best"]["performance"] * 200 # rough heuristic of when to quit
        parameter_population = fp_sample(fitness,
                                         parameter_population,
                                         sample_runs,
                                         performance_limit,
                                         num_threads)
        parameter_population = fp_cull(parameter_population, fitness_cutoff)
    if outpath:
        with open(outpath, 'w') as outfile:
            json.dump(tracking, outfile)
    else:
        pprint.pprint(tracking)
        
def main():
    """Run the mainloop and track fitness between generations."""
    parser = ArgumentParser()
    parser.add_argument("genetic_string")
    parser.add_argument("--meta", action="store_true", help="Enable the find-parameters meta-algorithm.")
    parser.add_argument("--meta-gen", type=int, default=10, dest="gens",
                        help="How many generations to run find-parameters for.")
    parser.add_argument("--meta-samples", type=int, default=5, dest="samples",
                        help="How many sample runs to do in find-parameters.")
    parser.add_argument("--meta-population", type=int, default=5, dest="mpop",
                        help="What population of parameters to use in find-parameters.")
    parser.add_argument("--meta-fitness", type=float, default=90, dest="mfp",
                        help="What fitness percentile from 1 to 99.999 to cull per generation of find-parameters.")
    parser.add_argument("--meta-threads", type=int, default=1, dest="threads",
                        help="How many threads to run find-parameters with.")
    parser.add_argument("--meta-output", type=str, default=None, dest="outpath",
                        help="Optional filepath to write history and results file to.")
    parser.add_argument("--fitness-percentile", type=float, default=90, dest="fp",
                        help=("Number from 1 to 99.999... specifying the " 
                              "fitness cutoff per generation."))
    parser.add_argument("--population", type=int, default=100, dest="pop",
                        help=("An integer number representing the population size to use."))
    parser.add_argument("--initial-ss", type=int, default=50, dest="string_size",
                        help="Integer representing initial strength length.")
    parser.add_argument("--test-performance", action="store_true",
                        dest="test_perf", help="Run cprofile on genefind")
    arguments = parser.parse_args()
    if arguments.meta:
        find_parameters(arguments.genetic_string,
                        generations=arguments.gens,
                        sample_runs=arguments.samples,
                        population=arguments.mpop,
                        fitness_percentile=arguments.mfp,
                        num_threads=arguments.threads,
                        outpath=arguments.outpath)
    else:                
        run_genefind(arguments.genetic_string,
                     arguments.fp,
                     arguments.pop,
                     arguments.string_size)
    
    
# Should really add an option to test performance
if __name__ == '__main__':
    main()
    
    
    
