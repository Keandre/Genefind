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
import random
import cProfile

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

#(???)
actg = "ACTG"

def random_delete_n(stringIn, n=1, consecutive=False):
    """Randomly delete n characters from a string."""
    stringOut = stringIn
    if len(stringOut)<n+1: return stringOut
    if not consecutive:
        for i in range(n):
            counter = random.randrange(len(stringOut))
            stringOut = stringOut[:counter]+stringOut[counter+1:]
        return stringOut
    else:
        counter = random.randrange((len(stringOut)-n+1))
        stringOut = stringOut[:counter]+stringOut[counter+n:]
        return stringOut

def random_add_n(stringIn, n=1, consecutive=False):
    """Randomly add n characters to a string."""
    stringOut = stringIn
    if len(stringOut)<n+1: return stringOut
    if not consecutive:
        for i in range(n):
            counter = random.randrange((len(stringOut)))
            stringOut = stringOut[:counter]+random.choice(actg)\
                        +stringOut[counter:]
        return stringOut
    else:
        counter = random.randrange((len(stringOut)-n+1))
        for i in range(n):
            stringOut = stringOut[:counter]+random.choice(actg)\
                        +stringOut[counter:]
        return stringOut

def RS(stringIn):
    stringOut = stringIn
    if len(stringOut)<2: return stringOut
    ct1 = random.randrange((len(stringOut)))
    ct2 = random.randrange((len(stringOut)))
    if ct1>ct2: ct2, ct1 = ct1, ct2
    stringOut = stringOut[:ct1]\
                +stringOut[ct2:ct2+1]\
                +stringOut[ct1+1:ct2]\
                +stringOut[ct1:ct1+1]\
                +stringOut[ct2+1:]
    return stringOut

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


# This might be confusing. Basically we're defining these operations in the global
# scope, and then executing them using eval(). This is about 1/3 slower than a native
# function call btw according to cProfiler.
operations = ["random_shift(string,3)","random_shift(string,1)","random_shift(string,-1)",
              "random_shift(string,-3)", "RS(string)",
              "random_add_n(string)","random_add_n(string,3)","random_add_n(string,3,True)",\
              "random_delete_n(string)","random_delete_n(string,3)","random_delete_n(string,3,True)"]

def Mod2Run(listIn, population):##pop gen
    """Generate the population. Here we randomly pull from a set of predefined
    operations to mutate our existing strings, then execute the operation."""
    #listIn = [a,b]
    mutCount = 0
    countProc = 0
    listOut = []
    listOut += listIn
    while len(listOut)<=population:
        string = listOut[countProc]
        op = random.choice(operations)
        listOut.append(eval(op))
        countProc+=1
    return listOut

def Mod3Run(string1, string2):##string comparison
    if len(string1)<2: return 0.0
    difference = sum(1 if x == y else 0 for x, y in \
                     zip(string1, string2)) \
                     * 100 / len(string2)
    lengthDif = min(len(string1)/len(string2),len(string2)/len(string1))
    difference *= lengthDif
    return difference

def Mod4Run(listIn, fitness, fitness_percentile):##pop cull
    listP = []
    for i in listIn:
        listP.append([Mod3Run(i,fitness),i])
    listP.sort()
    fitness_cutoff = -round(len(listP) * fitness_percentile)
    return listP[fitness_cutoff:]

def Mod1Run(listIn, trackIn={}):##iterator
    trackOut = trackIn
    bestE = listIn[0]
    listOut = []
    ##print(len(listIn))
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
              ,"String :\n", bestE[1])
        return listOut, trackOut
    return listOut, trackOut

def initial():
    output = ''
    for i in range(random.randint(20,50)):
        output+=random.choice(actg)
    return output

##base program
print("Modules loaded.")
print(\
    """
********************
"GeneFind: A Designspace Searcher"
********************
""")


def main():
    """Run the mainloop and track fitness between generations."""
    parser = ArgumentParser()
    parser.add_argument("genetic_string")
    parser.add_argument("--fitness-percentile", type=float, default=90, dest="fp",
                        help=("Number from 1 to 99.999... specifying the " 
                              "fitness cutoff per generation."))
    parser.add_argument("--population", type=int, default=100, dest="pop",
                        help=("An integer number representing the population size to use."))
    arguments = parser.parse_args()
    lis = [initial()]
    inputIsSafe = 0
    fitness = arguments.genetic_string
    fitness_cutoff = (100 + -arguments.fp) / 100
    track = {"generationCount":0, 'currentFit':0.0, \
             'lastFit':0.0, 'bestFit':0.0}
    while 1:
        lis, track = Mod1Run(lis, track)
        lis = Mod2Run(lis, arguments.pop)
        lis = Mod4Run(lis,fitness, fitness_cutoff)
        if track["bestFit"] == 100:
            break
        if track["generationCount"]%10 == 0 or track["bestFit"]>track["lastFit"]:
            # This should really be put under a -v option
            print("Generation ",track["generationCount"],"Fitness ", track["bestFit"])

# Should really add an option to test performance
main()
#cProfile.run('main()')
    
    
    
    
