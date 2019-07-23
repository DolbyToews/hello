import numpy
import matplotlib
import pylab
import random
import math
n=1000 #number of timesteps taken
k=2 #number of "type 1" particles spawning in area one
c=0 #number of "type 2" particles spawning in area one
g=0 #number of "type 1" particles spawning in area 2
f=2 #number of "type 2" particles spawning in area 2
d=2 #number of dimensions
E = numpy.zeros(k + c + g + f)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
a=0#lower bound for x and x-randomization
b=10#upper bound for x and x-randomization
h=7#separates areas into area one and area two, and sets y inter-area boundrary
y=4#lower bound for y and y-randomization
z=10#upper bound for y and y randomization
Trans=([[8, 7]])#need to test that this rejects stuff not in here, also, could replace "7" with "h"
T=1#thermal energy test value
DC1 = ([a], [b])#sets upper and lower boundaries at y=1 y= 5(divider) and y=10
#as well as boundaries at x=5x=10 x=1(unused due to shift of position)
DC3 = ([y], [h], [z])
A1x = numpy.random.randint(low=a + 1, high=b - 1, size =((k + c), 1))#randomization for area one x-bound
A1y = numpy.random.randint(low=y + 1, high=h - 1, size =((k + c), 1))#randomization for area one y-bound
A1 = numpy.concatenate((A1x, A1y), axis= 1)#matrix of positions of area one
A2x = numpy.random.randint(low=a + 1, high=b - 1, size =((g + f), 1))#randomization for area two x-bound
A2y = numpy.random.randint(low=h + 1, high=z - 1, size =((g + f), 1))#randomization for area two y-bound
A2 = numpy.concatenate((A2x, A2y), axis= 1)#matrix of positions of area one
if d == 1:#sets possible moves for each dimension
    M = numpy.array([(-1),(0),(1)])
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])
P = numpy.concatenate((A1, A2), axis= 0)#combines A1 and A2 ,as well as T1 and T2 particles,into a single matrix
#note that particles can't spawn on boundaries, but can spawn on eachother  
x = numpy.hsplit(P, d)
cool = x[0]
swell = x[1]
 #   alright = x[2]
s = numpy.vsplit(x[1], [k, (k + c), (k + c + g), (k + c + g + f)])
daf = s[0]
paf = s[1]
qaf = s[2]
zaf = s[3]
arealow = ((b - 1) - (a + 1)) * ((h) - (y + 1))
areahigh = ((b - 1) - (a + 1)) * ((z - 1) - (h))
fileout = open ("coordinatesold.txt", "w")
for i in range(0, n):
    type1low = 0
    type1high = 0
    type2low = 0
    type2high = 0
    j = random.randint(0, (k + c + g + f) - 1)#this picks which particle moves
    P_old = P[j].copy() #saves old position in case the new move is rejected
    P[j] = P[j] + random.choice(M) #the new move
    E_old = E[j] #saves old energy in case the new move is rejected
    for q in range(0, (k + c + g + f)): #includes all particles of all types
        if q == j: #these two lines make sure that the moved particle doesn't check against itself
            continue
        DC2 = P[q]
        if numpy.array_equal(P[j], DC2):#note, this whole sequence will not work if there's only one particle
            #it breaks because it doesn't go to the second step and thus doesn't block the boundaries
            o = 1 #this may or may not be a feature incorporated in the future
            E[j] = math.inf
            break 
#if the particle moves into the same position as another particle, it's position is infinite and the loop stops there
        elif Trans in cool[j] and h in swell[j]:
            E[j] = 0
            break
        elif DC1 in cool[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
            E[j] = math.inf
            break
        elif DC3 in swell[j]:#this isn't the best solution, as it runs through everything a couple times, however, it does have 100% accuracy
            E[j] = math.inf
            break
        #elif d == 3:
         #   if DC1 in alright[j]:
          #      E[j] = math.inf
           #     break
        else:
            E[j] = 0
    if E[j] > E_old:
        val = numpy.random.uniform(low=0, high=1, size=1)
        if val > (math.exp((E_old - E[j])/(T))):#equation will need to be changed later(rn inaccurate)
            #T represents thermal energy, unsure how we want that represented(Kb * T)
            #overall this will probably reject most moves
            P[j] = P_old
            E[j] = E_old
        #if the energy's lower we automatically accept it
        type1 = numpy.concatenate((daf, qaf), axis= 0)
        type2 = numpy.concatenate((paf, zaf), axis= 0)
    for q in range(0, (k + g)): #prints out each particle into the output file
        if type1[q] < h:
            type1low = type1low + 1
        elif type1[q] >= h:#maybe make it so that transporters aren't included in the area
            type1high = type1high + 1
    for q in range(0, (c + f)):
        if type2[q] < h:
            type2low = type2low + 1
        elif type2[q] >= h:
            type2high = type2high + 1
    for q in range(0, (k + c + g + f)): #prints out each particle into the output file
        print(*P[q], end = " ", file= fileout)
        if q == ((k + c + g + f) - 1): #checks if we're at the last particle
            print(file= fileout)
    densitylow = type1low + type2low 
    densityhigh = type1high + type2high
    print("type1low=", type1low/arealow, "type2low=", type2low/arealow, "lowdensity=", densitylow/arealow, file= fileout)
    print("type1high=", type1high/areahigh, "type2high=", type2high/areahigh, "highdensity=", densityhigh/areahigh, file= fileout)
fileout.close()   
