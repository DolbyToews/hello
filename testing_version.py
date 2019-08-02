import numpy
import matplotlib
import pylab
import random
import math
n=185 #number of timesteps taken
k=2 #number of "type 1" particles spawning in area one
c=0 #number of "type 2" particles spawning in area one
g=0 #number of "type 1" particles spawning in area 2
f=2 #number of "type 2" particles spawning in area 2
d=2 #number of dimensions
E = numpy.zeros(k + c + g + f)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
a=6#lower bound for x and x-randomization
b=13#upper bound for x and x-randomization
h=7#separates areas into area one and area two, and sets y inter-area boundrary
uph=h+1
y=4#lower bound for y and y-randomization
z=11#upper bound for y and y randomization
Trans=([7], [8])#x-coords on the cell membrane which act as holes
Special_Trans1=((9), (11))
Special_Trans2=((10), (12))
T=1#thermal energy test value
DC1 = ([a], [b])#sets upper and lower boundaries at x=1  and x=10
DC3 = ([y], [h], [uph], [z])#sets upper and lower y boundaries and the cell membrane level
DC4 = ([h], [uph])
A1x = numpy.random.randint(low=a + 1, high=b, size =((k + c), 1))#randomization for area one x-bound
A1y = numpy.random.randint(low=y + 1, high=h, size =((k + c), 1))#randomization for area one y-bound
A1 = numpy.concatenate((A1x, A1y), axis= 1)#matrix of positions of area one
A2x = numpy.random.randint(low=a + 1, high=b, size =((g + f), 1))#randomization for area two x-bound
A2y = numpy.random.randint(low=uph + 1, high=z, size =((g + f), 1))#randomization for area two y-bound
A2 = numpy.concatenate((A2x, A2y), axis= 1)#matrix of positions of area two
if d == 1:#sets possible moves for each dimension
    M = numpy.array([(-1),(0),(1)])
if d == 2:
    M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
if d == 3:
    M = numpy.array([(-1,-1,-1),(-1,-1,0),(-1,-1,1),(-1,0,-1),(-1,0,0),(-1,0,1),(-1,1,-1),(-1,1,0),(-1,1,1),(0,-1,-1),(0,-1,0),(0,-1,1),(0,0,-1),(0,0,0),(0,0,1),(0,1,-1),(0,1,0),(0,1,1),(1,-1,-1),(1,-1,0),(1,-1,1),(1,0,-1),(1,0,0),(1,0,1),(1,1,-1),(1,1,0),(1,1,1)])
P = numpy.concatenate((A1, A2), axis= 0)#combines A1 and A2 ,as well as T1 and T2 particles,into a single matrix
#note that particles can't spawn on boundaries, but can spawn on eachother  
x = numpy.hsplit(P, d)
x_cord = x[0]
y_cord = x[1]
arealow = ((b - 1) - (a + 1)) * ((h) - (y + 1))#calculates the total area of the lower compartment
areahigh = ((b - 1) - (a + 1)) * ((z - 1) - (h))#calculates the total area of the higher compartment
Ising = numpy.zeros((len(Trans), 2))
Idontsing = numpy.zeros((len(Special_Trans1), 3))
for u in range(0, len(Trans)):
    Ising[u] = numpy.zeros((1, 2))
for u in range(0, len(Special_Trans1)):
    Idontsing[u] = numpy.zeros((1, 3))
upspin = numpy.ones(2)
downspin = -1 * numpy.ones(2)
upone = numpy.ones(1)
downone = -1 * numpy.ones(1)
uponedownone = numpy.concatenate((upone, downone), axis= None)
downoneupone = numpy.concatenate((downone, upone), axis= None)
fileout = open ("coordinatesold.txt", "w")
print(P)
for i in range(0, n):
    type1low = 0
    type1high = 0
    type2low = 0
    type2high = 0
    y_cord_old = x[1].copy()
    j = random.randint(0, (k + c + g + f) - 1)#this picks which particle moves
    P_old = P[j].copy() #saves old position in case the new move is rejected
    P[j] = P[j] + random.choice(M) #the new move
    E_old = E[j] #saves old energy in case the new move is rejected
    I_old = Idontsing.copy()
    for q in range(0, (k + c + g + f)): #includes all particles of all types
        if q == j: #these two lines make sure that the moved particle doesn't check against itself
            continue
        DC2 = P[q]#makes sure that particles dont collide
        if numpy.array_equal(P[j], DC2):#note, this whole sequence will not work if there's only one particle
            #it breaks because it doesn't go to the second step and thus doesn't block the boundaries
            #o = 1 #this may or may not be a feature incorporated in the future
            E[j] = math.inf#energy set to infinite if they collide
            break #if the particle moves into the same position as another particle, it's energy is infinite and the loop stops there
        elif Trans in x_cord[j] and DC4 in y_cord[j]:#tells us if particle is on a transporter
            u = Trans.index(x_cord[j])
            if uph in y_cord_old[j] and h in y_cord[j]:
                if upspin in Ising[u]: #rejects particles heading "down" if the transporter is receiving particles moving up
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j)
                Ising[u] = upspin #sets transporter to "up" if this is the first particle to pass through
            if h in y_cord_old[j] and uph in y_cord[j]:
                if downspin in Ising[u]: #rejects particle heading "up" if the transporter is receiving particles moving "down"
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j)
                Ising[u] = downspin#ses transporter to "down" if this is the first particle to pass through
            break
        elif Special_Trans1 in x_cord[j] and DC4 in y_cord[j]:#sets which"end" of transporter is altered
            u = Special_Trans1.index(x_cord[j])#sets which transporter is altered
            funandgames = numpy.split(Idontsing[u], [2])#sets up so that we can alter the specific transporter end
            fun_old = funandgames[1].copy()#might remove this at some point, not sure if it's needed anymore or not
            if uph in y_cord_old[j] and h in y_cord[j]:
                if numpy.array_equal(upspin, funandgames[0]):
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j)
                elif numpy.array_equal(uponedownone, funandgames[0]):
                    E[j] = 0
                    funandgames[0] = downoneupone
                    Idontsing[u] = numpy.concatenate((funandgames[0], fun_old), axis= None)
                    print("glitch", funandgames[0], fun_old)
                elif numpy.array_equal(downoneupone, funandgames[0]):
                    E[j] = 0
                    funandgames[0] = uponedownone
                    Idontsing[u] = numpy.concatenate((funandgames[0], fun_old), axis= None)
                    print("glitch", funandgames[0], fun_old)
                elif numpy.array_equal(downspin, funandgames[0]):
                    E[j] = 0
                    funandgames[0] = upspin
                    Idontsing[u] = numpy.concatenate((funandgames[0], fun_old), axis= None)
                    print("glitch", funandgames[0], fun_old)
                else:
                    E[j] = 0
                    funandgames[0] = upspin
                    Idontsing[u] = numpy.concatenate((funandgames[0], fun_old), axis= None)
            if h in y_cord_old[j] and uph in y_cord[j]:
                if numpy.array_equal(downspin, funandgames[0]):
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j)
                elif numpy.array_equal(uponedownone, funandgames[0]):
                    E[j] = 0    
                    funandgames[0] = downoneupone
                    Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
                    print("glitch", funandgames[0], funandgames[1])
                elif numpy.array_equal(downoneupone, funandgames[0]):
                    E[j] = 0
                    funandgames[0] = uponedownone
                    Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
                    print("glitch", funandgames[0], funandgames[1])
                elif numpy.array_equal(upspin, funandgames[0]):
                    E[j] = 0
                    funandgames[0] = downspin
                    Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
                    print("glitch", funandgames[0], funandgames[1])
                else:
                    E[j] = 0
                    funandgames[0] = downspin
                    Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None) 
            break
        elif Special_Trans2 in x_cord[j] and DC4 in y_cord[j]:
            u = Special_Trans2.index(x_cord[j])
            funny = numpy.split(Idontsing[u], [1])
            funny_old = funny[0].copy()#might remove this at some point, not sure if it's needed anymore or not
            if uph in y_cord_old[j] and h in y_cord[j]:
                if numpy.array_equal(upspin, funny[1]):
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j) 
                elif numpy.array_equal(uponedownone, funny[1]):
                    E[j] = 0
                    funny[1] = downoneupone
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                elif numpy.array_equal(downoneupone, funny[1]):
                    E[j] = 0
                    funny[1] = uponedownone
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                elif numpy.array_equal(downspin, funny[1]):
                    E[j] = 0
                    funny[1] = upspin
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                else:
                    E[j] = 0
                    funny[1] = upspin
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
            if h in y_cord_old[j] and uph in y_cord[j]:
                if numpy.array_equal(downspin, funny[1]):
                    E[j] = math.inf
                    print("transporter closed")
                    print("rejected=", P[j], j)
                elif numpy.array_equal(uponedownone, funny[1]):
                    E[j] = 0    
                    funny[1] = downoneupone
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                elif numpy.array_equal(downoneupone, funny[1]):
                    E[j] = 0
                    funny[1] = uponedownone
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                elif numpy.array_equal(upspin, funny[1]):
                    E[j] = 0
                    funny[1] = downspin
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
                    print("glitch", funny_old, funny[1])
                else:
                    E[j] = 0
                    funny[1] = downspin
                    Idontsing[u] = numpy.concatenate((funny_old, funny[1]), axis= None)
            break
        elif DC1 in x_cord[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
            E[j] = math.inf#energy set to infinite if particle runs into x boundaries
        elif DC3 in y_cord[j]:#this isn't the best solution, as it runs through everything a couple times, however, it does have 100% accuracy
            E[j] = math.inf#energy set to infinite if particle runs into y boundaries and membrane
        else:
            E[j] = 0
    if E[j] > E_old:#if the energy's lower or equal we automatically accept it
        val = numpy.random.uniform(low=0, high=1, size=1)
        if val > (math.exp((E_old - E[j])/(T))):#equation will need to be changed later(rn inaccurate)
            #T represents thermal energy, unsure how we want that represented(Kb * T)
            P[j] = P_old
            E[j] = E_old
            Idontsing = I_old
    for q in range(0, (k + c + g + f)): #prints out each particle into the output file
        print(*P[q], end = " ", file= fileout)
        if q == ((k + c + g + f) - 1): #checks if we're at the last particle
            print(file= fileout)
    s = numpy.vsplit(x[1], [k, (k + c), (k + c + g), (k + c + g + f)])#splits up y-coordinates of particles
    A1T1 = s[0]#y-coord of type 1 particles spawning in the lower area
    A1T2 = s[1]#y-coord of type 2 particles spawning in the lower area
    A2T1 = s[2]#y-coord of type 1 particles spawning in the higher area
    A2T2 = s[3]#y-coord of type 2 particles spawning in the higher area
    type1 = numpy.concatenate((A1T1, A2T1), axis= 0) #aggregates all type one particles regardless of spawn
    type2 = numpy.concatenate((A1T2, A2T2), axis= 0) #aggregates all type two particles regardless of spawn
    for q in range(0, (k + g)): 
        if type1[q] < h: #tells us how many type one particles are in area one
            type1low = type1low + 1
        elif type1[q] >= h:#tells us how many type one particles are in area two
            #maybe make it so that transporters aren't included in area two, not sure yet
            type1high = type1high + 1
    for q in range(0, (c + f)):
        if type2[q] < h:#tells us how many type two particles are in area one
            type2low = type2low + 1
        elif type2[q] >= h: #tells us how many type two particles are in area two
            type2high = type2high + 1
    densitylow = type1low + type2low #number of particles in the lower compartment
    densityhigh = type1high + type2high #number of particles in the higher compartment
    print("type1low=", type1low/arealow, "type2low=", type2low/arealow, "lowdensity=", densitylow/arealow, file= fileout)
    print("type1high=", type1high/areahigh, "type2high=", type2high/areahigh, "highdensity=", densityhigh/areahigh, file= fileout)
    print(j, P[j], Idontsing)
fileout.close() 