import numpy
import matplotlib
import pylab
import random

n=10 #number of timesteps taken
k=3 #number of particles
d=2 #number of dimensions

P = numpy.random.randint(low=0, high=10, size =(k, d))
if d == 1:
    M = numpy.array([(-1),(0),(1)])
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])
   
fileout = open ("coordinatesold.txt", "w")
for i in range(0, n):
    j = random.randint(0, k - 1)#this picks which particle moves
    #print(i, P[j])
    P[j] = P[j] + random.choice(M)
    #for j in range(0, k):
    #if val == j:
     #   P[j] = P[j] + random.choice(M)
    #print(i, P[j], "particle=", j + 1)
    for q in range(0, k):
        print(*P[q], end = " ", file= fileout)
        if q == k - 1: #checks if we're at the last particle
            print(file= fileout)
            #change this to write an arbitrary number of particles as columns next to eachother(how?)
            #Same for MDS
            #Make Diffusion Constant whatever the MDS is divided by 2d
            #selects the "x" and "y" coordinates of each P[j] for calculation of MSD
fileout.close()
