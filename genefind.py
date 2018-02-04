#mutations:
#dshift by [val={1,3}]
#drandom delete 1
#drandom delete 3 consecutive
#drandom delete 3 random (essentially RD1 three times)
#drandom add 1
#drandom add 3 consecutive
#drandom add 3 random
#drandom switch

import random

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

def RDV(stringIn, val=1, consecutive=False):
    stringOut = stringIn
    if len(stringOut)<val+1: return stringOut
    if not consecutive:
        for i in range(val):
            counter = random.choice(range(len(stringOut)))
            stringOut = stringOut[:counter]+stringOut[counter+1:]
        return stringOut
    else:
        counter = random.choice(range(len(stringOut)-val+1))
        stringOut = stringOut[:counter]+stringOut[counter+val:]
        return stringOut

def RAV(stringIn, val=1, consecutive=False):
    stringOut = stringIn
    if len(stringOut)<val+1: return stringOut
    if not consecutive:
        for i in range(val):
            counter = random.choice(range(len(stringOut))) # Inefficient
            stringOut = stringOut[:counter]+random.choice(actg)\
                        +stringOut[counter:]
            #print(stringOut)
        return stringOut
    else:
        counter = random.choice(range(len(stringOut)-val+1))
        for i in range(val):
            stringOut = stringOut[:counter]+random.choice(actg)\
                        +stringOut[counter:]
            #print(stringOut)
        return stringOut

def RS(stringIn):
    stringOut = stringIn
    if len(stringOut)<2: return stringOut
    ct1 = random.choice(range(len(stringOut)))
    ct2 = random.choice(range(len(stringOut)))
    if ct1>ct2: ct2, ct1 = ct1, ct2
    stringOut = stringOut[:ct1]\
                +stringOut[ct2:ct2+1]\
                +stringOut[ct1+1:ct2]\
                +stringOut[ct1:ct1+1]\
                +stringOut[ct2+1:]
    return stringOut

def RSH(stringIn, val):
    if len(stringIn)<val+1: return stringIn
    stringOut = stringIn[val:]+stringIn[:val]
    return stringOut

operations = ["RSH(string,3)","RSH(string,1)","RSH(string,-1)","RSH(string,-3)",\
              "RS(string)","RAV(string)","RAV(string,3)","RAV(string,3,True)",\
              "RDV(string)","RDV(string,3)","RDV(string,3,True)"]

def Mod2Run(listIn):##pop gen
    #listIn = [a,b]
    mutCount = 0
    countProc = 0
    listOut = []
    listOut += listIn
    while len(listOut)<=100:
        ##print(countProc, len(listOut))
        string = listOut[countProc]
        op = random.choice(operations)
        listOut.append(eval(op))
        countProc+=1
        #print(countProc, listOut[len(listOut)-1], op, len(listOut[len(listOut)-1]))
    #for i in listOut:print(countProc, i)
    return listOut

def Mod3Run(string1, string2):##string comparison
    if len(string1)<2: return 0.0
    difference = sum(1 if x == y else 0 for x, y in \
                     zip(string1, string2)) \
                     * 100 / len(string2)
    lengthDif = min(len(string1)/len(string2),len(string2)/len(string1))
    difference *= lengthDif
    return difference

def Mod4Run(listIn, fitness):##pop cull
    listP = []
    for i in listIn:
        listP.append([Mod3Run(i,fitness),i])
    listP.sort()
    listP.reverse()
    while len(listP)>10:
        listP.pop()
    listOut = listP
##    for i in listP:
##        listOut.append(i[1])
    return listOut

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
lis = [initial()]
inputIsSafe = 0
##while inputIsSafe:
    ##fitness = input("Input comparison string to be found: ")
fitness = input("Input comparison string to be found: ")
track = {"generationCount":0, 'currentFit':0.0, \
                    'lastFit':0.0, 'bestFit':0.0}
while 1:
    lis, track = Mod1Run(lis, track)
##    if track["generationCount"]<10 or track["generationCount"]%100 == 0:
##        print(track["bestFit"], lis)
    lis = Mod2Run(lis)
##    if track["generationCount"]<10 or track["generationCount"]%100 == 0:
##        print(track["bestFit"], lis)
    lis = Mod4Run(lis,fitness)
    if track["bestFit"] == 100:
        break
    if track["generationCount"]%10 == 0 or track["bestFit"]>track["lastFit"]:
        print("Generation ",track["generationCount"],"Fitness ", track["bestFit"])
        ##print("*****")
    
    
    
    
