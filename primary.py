import numpy
import matplotlib
import pylab
import random
n=10 #number of timesteps taken
k=3 #number of particles
d=2 #number of dimensions
E = numpy.zeros(n)#default energy for now, may change to a equation and move placement of it
a=-40#lower bound for randomization, and overall boundaries
b=40#upper bound for randomization, and overall boundaries
DC1 = ([a], [b])
P = numpy.random.randint(low=a + 1, high=b - 1, size =(k, d))#matrix of positions
if d == 1:
    M = numpy.array([(-1),(0),(1)])
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])  
fileout = open ("coordinatesold.txt", "w")
for i in range(0, n):
    j = random.randint(0, k - 1)#this picks which particle moves
    P[j] = P[j] + random.choice(M)
    if d in range(1, 4):#we shouldn't need to need "d == 1", lines 25 through 31 should work fine, but it doesn't work fine
        if DC1 in P[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
            E[i] = 1     
            #print("e=", E[i])#change output if we want to export to text file
        elif DC1 not in P[j]:
            E[i] = 0
    if E[i] > E[i - 1]:
        P[j][i] = P[j][i - 1]#problem is this
    elif E[i] <= E[i - 1]:
        P[j][i] = P[j][i]# and also this
    for q in range(0, k):
        print(*P[q], end = " ", file= fileout)
        if q == k - 1: #checks if we're at the last particle
            print(file= fileout)
            #change this to write an arbitrary number of particles as columns next to eachother(how?)
fileout.close()
#D = numpy.genfromtxt("coordinatesold.txt")
#for q in range(1, (2 * k), 2): #this makes graphs of the path of each particle
 #   pylab.title("Random Walk ($n = " + str(n) + "$ steps)") 
  #  pylab.plot(D[:, q - 1], D[:, q]) 
   # pylab.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600) 
    #pylab.show() 
