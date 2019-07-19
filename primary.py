import numpy
import matplotlib
import pylab
import random
import math
n=10 #number of timesteps taken
k=2 #number of "type 1" particles spawning in area one
c=0 #number of "type 2" particles spawning in area one
g=0 #number of "type 1" particles spawning in area 2
f=2 #number of "type 2" particles spawning in area 2
d=2 #number of dimensions
E = numpy.zeros(k + c + g + f)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
a=0#lower bound for overall boundaries and second square randomization
b=10#upper bound for overall boundaries
h=5#separates areas into area one and area two, and sets x=5 boundrary
T=1#thermal energy test value
DC1 = ([a], [h], [b])#sets upper and lower boundaries at y=1 y= 5(divider) and y=10
#as well as boundaries at x=5x=10 x=1(unused due to shift of position)
A = numpy.random.randint(low=a + 1, high=h - 1, size =((k + c), d))#matrix of positions of area one
A2 = numpy.random.randint(low=h + 1, high=b - 1, size =((g + f), d))#matrix of positions of area two
A5 = ((h + 1) - (a + 1)) * numpy.ones((k + c, 1)) #shifting every particle in area 1 along the x-axis
A0 = numpy.zeros((k + c, 1)) #keeps shift from affecting the y-axis
if d == 1:#sets possible moves for each dimension
    M = numpy.array([(-1),(0),(1)])
    A1 = A
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
    A1 = A + numpy.concatenate((A5, A0), axis= 1) #finalizes move on x-axis
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])
    A1 = A + numpy.concatenate((A5, A5, A0), axis= 1) #moves on x and y axes when in 3D
P = numpy.concatenate((A1, A2), axis= 0)#combines A1 and A2 ,as well as T1 and T2 particles,into a single matrix
#note that particles can't spawn on boundaries, but can spawn on eachother  
fileout = open ("coordinatesold.txt", "w")
for i in range(0, n):
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
        else:
            o = 0
            if DC1 in P[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
                E[j] = math.inf
            elif DC1 not in P[j]:#this isn't the best solution, as it runs through everything a couple times, however, it does have 100% accuracy
                E[j] = 0         
    if E[j] > E_old:
        val = numpy.random.uniform(low=0, high=1, size=1)
        if val > (math.exp((E_old - E[j])/(T))):#equation will need to be changed later(rn inaccurate)
            #T represents thermal energy, unsure how we want that represented(Kb * T)
            #overall this will probably reject most moves
            P[j] = P_old
            E[j] = E_old
        #if the energy's lower we automatically accept it
    for q in range(0, (k + c + g + f)): #prints out each particle into the output file
        print(*P[q], end = " ", file= fileout)
        if q ==  ((k + c + g + f) - 1): #checks if we're at the last particle
            print(file= fileout)
fileout.close()    
