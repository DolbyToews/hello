import numpy
import matplotlib
import pylab
import random
import math
n=1000#number of timesteps taken
k=2 #number of "type 1" particles spawning in area one
c=0 #number of "type 2" particles spawning in area one
g=0 #number of "type 1" particles spawning in area 2
f=2 #number of "type 2" particles spawning in area 2
Reject_Type_One = range(0, k) or range(k + c, k + c + g)
Reject_Type_Two = range(k, k + c) or range(k + c + g, k + c + g + f)
a=9#lower bound for x and x-randomization
b=13#upper bound for x and x-randomization
h=7#separates areas into area one and area two, and sets y inter-area boundrary
uph=h+1
y=5#lower bound for y and y-randomization
z=10#upper bound for y and y randomization
Trans=()#x-coords on the cell membrane which act as holes
Special_Trans1=((9), (11))#sets part of transporter which rejects type one particles
Special_Trans2=((10), (12))#sets part of transporter which rejects type two particles
#in order for this to work, they must have the same number of entries
T=1#thermal energy test value
DC1 = ([a], [b])#sets upper and lower boundaries at x=1  and x=10
DC3 = ([y], [h], [uph], [z])#sets upper and lower y boundaries and the cell membrane level
DC4 = ([h], [uph])
M = numpy.array([(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1),(1,0),(0,1),(1,1)])
A1x = numpy.random.randint(low=a + 1, high=b, size =((k + c), 1))#randomization for area one x-bound
A1y = numpy.random.randint(low=y + 1, high=h, size =((k + c), 1))#randomization for area one y-bound
A1 = numpy.concatenate((A1x, A1y), axis= 1)#matrix of positions of area one
A2x = numpy.random.randint(low=a + 1, high=b, size =((g + f), 1))#randomization for area two x-bound
A2y = numpy.random.randint(low=uph + 1, high=z, size =((g + f), 1))#randomization for area two y-bound
A2 = numpy.concatenate((A2x, A2y), axis= 1)#matrix of positions of area two
P = numpy.concatenate((A1, A2), axis= 0)#sets matrix of positions, with particles spawning in area two being below those spawning in area one
x = numpy.hsplit(P, 2)
x_cord = x[0]#sets matrix of all x coordinates of all particles
y_cord = x[1]#sets matrix of all y coordinates of all particles
E_young = 0#default energy for zero particles in a transporter
arealow = ((b - 1) - (a + 1)) * ((h) - (y + 1))#calculates the total area of the lower compartment
areahigh = ((b - 1) - (a + 1)) * ((z - 1) - (h))#calculates the total area of the higher compartment
Ising = numpy.zeros((len(Trans), 2))#matrix of spins for transporters acting as protein channels
Idontsing = numpy.zeros((len(Special_Trans1), 2))#matrix of spins for transporters acting as coupled symporters/antiporters
Transrand = len(Special_Trans1)#sets number of transporters acting as such
E = numpy.zeros((2 * Transrand) + k + c + g + f)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
o = 0
for u in range(0, len(Trans)):
    Ising[u] = numpy.zeros((1, 2))#sets the initial spin of each individual protein channel
upspin = numpy.ones(2)
downspin = -1 * numpy.ones(2)
upone = numpy.ones(1)
downone = -1 * numpy.ones(1)
uponedownone = numpy.concatenate((upone, downone), axis= None)
downoneupone = numpy.concatenate((downone, upone), axis= None)
Special = (Special_Trans1, Special_Trans2)
Sppecial = numpy.zeros((len(Special_Trans1), 2))
for u in range(0, len(Special_Trans1)):
    Sppecial[u] = [item[u] for item in Special]
for u in range(0, len(Special_Trans1)):
    Idontsing[u] = uponedownone
#sets starting transporter positions to be "open"
KbT=1
epsilion=1
negepsilion=-1
partcount=1
def func1():
    if DC1 in x_cord[j]:#keeps particles from leaving the system's boundaries
        E[j] = math.inf#energy set to infinite if particle runs into x boundaries
        print("bound", E[j], P[j])
    elif Trans in x_cord[j] and DC4 in y_cord[j]:#tells us if particle is on a transporter
        u = Trans.index(x_cord[j])
        if uph in y_cord_old[j] and h in y_cord[j]:
            if upspin in Ising[u]: #rejects particles heading "down" if the transporter is receiving particles moving up
                E[j] = math.inf#keeps particles out if protein channel is closed
                print("transporter closed")
                print("rejected=", P[j], j)
            Ising[u] = upspin #sets transporter to "up" if this is the first particle to pass through
        if h in y_cord_old[j] and uph in y_cord[j]:
            if downspin in Ising[u]: #rejects particle heading "up" if the transporter is receiving particles moving "down"
                E[j] = math.inf
                print("transporter closed")
                print("rejected=", P[j], j)
            Ising[u] = downspin#ses transporter to "down" if this is the first particle to pass through
    elif Special_Trans1 in x_cord[j] and DC4 in y_cord[j]:#sets which"end" of transporter is altered
        u = Special_Trans1.index(x_cord[j])#sets which transporter is altered
        funandgames = numpy.split(Idontsing[u], 2)#sets up so that we can alter the specific transporter end
        if h in y_cord[j]:
            if j in Reject_Type_One:
                E[j] = math.inf#rejects type one particles from entering the special transporters from this way
                print("type one rejected")
                print("rejected=", P[j], j)
            elif numpy.array_equal(upone, funandgames[1]):
                E[j] = negepsilion#if the spin is up, the particle has a negative energy
                print("transporter closed")
                print("rejected=", P[j], j)
            elif numpy.array_equal(downone, funandgames[1]):
                E[j] = 0
                print("glitch", funandgames[0], funandgames[1])
            else:
                E[j] = 0
        elif uph in y_cord[j]:
            if j in Reject_Type_One:
                E[j] = math.inf
                print("type one rejected")
                print("rejected=", P[j], j)
            elif numpy.array_equal(upone, funandgames[0]):
                E[j] = negepsilion
                print("transporter closed")
                print("rejected=", P[j], j)
            elif numpy.array_equal(downone, funandgames[0]):
                E[j] = 0
                print("glitch", funandgames[0], funandgames[1])
            else:
                E[j] = 0
    elif Special_Trans2 in x_cord[j] and DC4 in y_cord[j]:
        u = Special_Trans2.index(x_cord[j])
        funandgames = numpy.split(Idontsing[u], 2)
        if h in y_cord[j]:
            if j in Reject_Type_Two:
                E[j] = math.inf#rejects type two particles
                print("type two rejected")
                print("rejected=", P[j], j)
            elif numpy.array_equal(upone, funandgames[1]):
                E[j] = negepsilion
                print("transporter closed")
                print("rejected=", P[j], j)
            elif numpy.array_equal(downone, funandgames[1]):
                E[j] = 0
                print("glitch", funandgames[0], funandgames[1])
            else:
                E[j] = 0
        elif uph in y_cord[j]:
            if j in Reject_Type_Two:
                E[j] = math.inf
                print("type two rejected")
                print("rejected=", P[j], j)
            elif numpy.array_equal(upone, funandgames[0]):
                E[j] = negepsilion
                print("rejected=", P[j], j)
            elif numpy.array_equal(downone, funandgames[0]):
                E[j] = 0
                print("glitch", funandgames[0], funandgames[1])
            else:
                E[j] = 0 
    elif DC3 in y_cord[j]:#this isn't the best solution, as it runs through everything a couple times, however, it does have 100% accuracy
        E[j] = math.inf#energy set to infinite if particle runs into y boundaries and membrane
        print("bound", E[j], P[j])
    else:
        E[j] = 0#default energy is 0
def func2():#sets energy of transporter(for the end of transporter that rejects type one) when it is randomized
    u = j - (k + c + g + f)
    funandgames = numpy.hsplit(Idontsing[u], 2)
    for u in range(0, len(Special_Trans1)):
        if funandgames[0] in upone:
            if funandgames[1] in upone:
                E[j] = epsilion
            if funandgames[1] in downone:
                E[j] = 0
        if funandgames[0] in downone:
            if funandgames[1] in upone:
                E[j] = 0
            if funandgames[1] in downone:
                E[j] = epsilion
def func3():#sets energy of transporter(for the end of transporter that rejects type two) when it is randomized
    u = j - (Transrand + k + c + g + f)
    funandgames = numpy.hsplit(Idontsing[u], 2)
    for u in range(0, len(Special_Trans1)):
        if funandgames[1] in upone:
            if funandgames[0] in upone:
                E[j] = epsilion
            if funandgames[0] in downone:
                E[j] = 0
        if funandgames[1] in downone:
            if funandgames[0] in upone:
                E[j] = 0
            if funandgames[0] in downone:
                E[j] = epsilion
fileout = open ("coordinatesold.txt", "a")
print(P)
for i in range(0, n):
    E_young = 0
    type1low = 0#counts what particles of what types are in what area
    type1high = 0
    type2low = 0
    type2high = 0
    y_cord_old = x[1].copy()#saves old y-coordinate, may not be neccessary anymore
    j = random.randint(0, ((2 * Transrand) + k + c + g + f) - 1)#this picks which particle moves
    E_old = E.copy() #saves old energy in case the new move is rejected
    Ej_old = E[j].copy()
    I_old = Idontsing.copy()
    if j in range(0, (k + c + g + f)):
        P_old = P[j].copy() #saves old position in case the new move is rejected
        P[j] = P[j] + random.choice(M) #the new move
        for q in range(0, (k + c + g + f)): #includes all particles of all types
            if q == j: #these two lines make sure that the moved particle doesn't check against itself
                continue
            DC2 = P[q]#makes sure that particles dont collide
            if numpy.array_equal(P[j], DC2):#note, this whole sequence will not work if there's only one particle
                #it breaks because it doesn't go to the second step and thus doesn't block the boundaries
                #o = 1 #this may or may not be a feature incorporated in the future
                E[j] = math.inf#energy set to infinite if they collide
                o = 1
                print("particle collision")
                break #if the particle moves into the same position as another particle, it's energy is infinite and the loop stops there
            else:
                o = 0
    if j in range(0, (k + c + g + f)):
        if o == 0:
            func1()#particle only checks if it runs into boundaries if it doesn't run into another particle
        else:
            E[j] = math.inf#if it ran into another particle, then it's rejected
    elif j in range((k + c + g + f), (Transrand + k + c + g + f)):
        print("yay")
        u = j - (k + c + g + f)#j=4 and j=5 are on two separate transporters that reject type one particles
        funandgames = numpy.hsplit(Idontsing[u], 2)
        if funandgames[0] in upone:#sets energy of a change in transporter state
            print("cool", u, j)
            funandgames[0] = downone#changes transporter's state
            if funandgames[1] in downone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount#accounts for change in particle energy within transporter when it's changing
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
                cool = E[j].copy()
            elif funandgames[1] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
                cool = E[j].copy()
        elif funandgames[0] in downone:
            print("cool", u, j)
            funandgames[0] = upone
            if funandgames[1] in downone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
                cool = E[j].copy()
            elif funandgames[1] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
                cool = E[j].copy()
        print("E[j]=", E[j])
        Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
        print("lol", i, j, P, Idontsing)
        for j in range(0, (k + c + g + f)):
            func1()
        for j in range((k + c + g + f), ((Transrand) + k + c + g + f)):
            func2()
        for j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):
            func3()
    elif j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):#same as above for the other end of the transporter (yes, this is neccessary)
        print("yay")
        u = j - (Transrand + k + c + g + f)
        funandgames = numpy.hsplit(Idontsing[u], 2)
        if funandgames[1] in upone:
            print("cool", u, j)
            funandgames[1] = downone
            if funandgames[0] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
                cool = E[j].copy()
            elif funandgames[0] in downone:      
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
                cool = E[j].copy()
        elif funandgames[1] in downone:
            print("cool", u, j) 
            funandgames[1] = upone
            if funandgames[0] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
                cool = E[j].copy()
            elif funandgames[0] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
                cool = E[j].copy()
        print("E[j]=", E[j])
        Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
        print("lol", i, j, P, Idontsing)
        for j in range(0, (k + c + g + f)):
            func1()
        for j in range((k + c + g + f), ((Transrand) + k + c + g + f)):
            func2()
        for j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):
            func3()
        #these three functions calculate the energy of all particles and transporters after a transporter shift occurs 
    print("norm", i, P, E, Idontsing)
    if j in range(0, (k + c + g + f)): #resets to old state if the move is rejected
        if E[j] > Ej_old:#if the energy's lower or equal we automatically accept it
            val = numpy.random.uniform(low=0, high=1, size=1)
            if val > (math.exp((Ej_old - E[j])/(KbT))):#equation will need to be changed later(rn inaccurate)
                #T represents thermal energy, unsure how we want that represented(Kb * T)
                P[j] = P_old
                E = E_old
                Idontsing = I_old
    if j in range((k + c + g + f), ((2 * Transrand) + (k + c + g + f))):
        if cool > Ej_old:
            val = numpy.random.uniform(low=0, high=1, size=1)
            if val > (math.exp((Ej_old - cool)/(KbT))):#equation will need to be changed later(rn inaccurate)
                print("fat=", cool)
                print(Ej_old - cool)
                E = E_old
                Idontsing = I_old
            print("get rekt", i, P, E[j], Idontsing)
    for q in range(0, (k + c + g + f)): #prints out each particle into the output file
        print(*P[q], end = " ", file= fileout)
        if q == ((k + c + g + f) - 1): #checks if we're at the last particle
            print(file= fileout)
    s = numpy.vsplit(x[1], [k, (k + c), (k + c + g), (k + c + g + f)])#splits up y-coordinates of particles
    #if j in range(0, (k + c + g + f)):
     #   print(i, j, P[j], Idontsing)
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
fileout.close() 
fileout = open ("p.txt", "w")
for j in range(0, (k + c + g +f)):
    print(*P[j], file= fileout)
fileout.close()#creates a new file from which a simulation can be continued
