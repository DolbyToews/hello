import numpy
import matplotlib
import pylab
import random
import math
n=10 #number of timesteps taken
k=4 #number of particles
d=2 #number of dimensions
E = numpy.zeros(k)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
a=0#lower bound for overall boundaries
b=3#upper bound for overall boundaries
T=1#thermal energy test value
DC1 = ([a], [b])
P = numpy.random.randint(low=a + 1, high=b - 1, size =(k, d))#matrix of positions
#makes sure that particles don't spawn on the boundary
if d == 1:#sets possible moves for each dimension
    M = numpy.array([(-1),(0),(1)])
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])  
fileout = open ("coordinatesold.txt", "w")
for i in range(0, n):
    j = random.randint(0, k - 1)#this picks which particle moves
    P_old = P[j].copy()
    P[j] = P[j] + random.choice(M)
    E_old = E[j]
    for q in range(0, k):
        if q == j:
            continue
        DC2 = P[q]
        if numpy.array_equal(P[j], DC2):
            o = 1 #this may or may not be a feature incorporated in the future
            E[j] = math.inf
            break
        else:
            o = 0
            if DC1 in P[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
                E[j] = math.inf
            elif DC1 not in P[j]:
                E[j] = 0         
    if E[j] > E_old:
        val = numpy.random.uniform(low=0, high=1, size=1)
        if val > (math.exp((E_old - E[j])/(T))):#equation will need to be changed later(rn inaccurate)
            #T represents thermal energy, unsure how we want that represented(Kb * T)
            P[j] = P_old
            E[j] = E_old
        #if the energy's lower we automatically accept it
    for q in range(0, k):
        print(*P[q], end = " ", file= fileout)
        if q == k - 1: #checks if we're at the last particle
            print(file= fileout)
            #change this to write an arbitrary number of particles as columns next to eachother(how?)
fileout.close()
