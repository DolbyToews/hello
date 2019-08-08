import numpy
import matplotlib
import pylab
import random
import math
n=500000#number of timesteps taken
k=40 #number of "type 1" particles spawning in area one
c=0 #number of "type 2" particles spawning in area one
g=0 #number of "type 1" particles spawning in area 2
f=40 #number of "type 2" particles spawning in area 2
Reject_Type_One = range(0, k) or range(k + c, k + c + g)
Reject_Type_Two = range(k, k + c) or range(k + c + g, k + c + g + f)
a=-1#lower bound for x and x-randomization
b=11#upper bound for x and x-randomization
h=10#separates areas into area one and area two, and sets y inter-area boundrary
uph=h+1
y=-1#lower bound for y and y-randomization
z=21#upper bound for y and y randomization
Trans=()#x-coords on the cell membrane which act as holes
Special_Trans1=((9), (11))
Special_Trans2=((10), (12))
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
P = numpy.concatenate((A1, A2), axis= 0)
x = numpy.hsplit(P, 2)
x_cord = x[0]
y_cord = x[1]
E_young = 0
arealow = ((b - 1) - (a + 1)) * ((h) - (y + 1))#calculates the total area of the lower compartment
areahigh = ((b - 1) - (a + 1)) * ((z - 1) - (h))#calculates the total area of the higher compartment
Ising = numpy.zeros((len(Trans), 2))
Idontsing = numpy.zeros((len(Special_Trans1), 2))
Transrand = len(Special_Trans1)
E = numpy.zeros((2 * Transrand) + k + c + g + f)#default energy for now, may change to a equation and move placement of it
#creates a energy for each particle for each position
for u in range(0, len(Trans)):
    Ising[u] = numpy.zeros((1, 2))
upspin = numpy.ones(2)
downspin = -1 * numpy.ones(2)
upone = numpy.ones(1)
downone = -1 * numpy.ones(1)
uponedownone = numpy.concatenate((upone, downone), axis= None)
downoneupone = numpy.concatenate((downone, upone), axis= None)
Special = (Special_Trans1, Special_Trans2)
Sppecial = numpy.zeros((len(Special_Trans1), 2))
type1highs = numpy.zeros(n)#remove these four lines when not doing diffussion testing
type2highs = numpy.zeros(n)
type1lows = numpy.zeros(n)
type2lows = numpy.zeros(n)
lol = numpy.arange(n)
for u in range(0, len(Special_Trans1)):
    Sppecial[u] = [item[u] for item in Special]
for u in range(0, 2):
    Sppecial[u] = [item[u] for item in Special]
print(Sppecial[1])
for u in range(0, len(Special_Trans1)):
    Idontsing[u] = uponedownone
KbT=1
epsilion=1
negepsilion=-1
partcount=1
def func1():
    if Trans in x_cord[j] and DC4 in y_cord[j]:#tells us if particle is on a transporter
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
    elif Special_Trans1 in x_cord[j] and DC4 in y_cord[j]:#sets which"end" of transporter is altered
        u = Special_Trans1.index(x_cord[j])#sets which transporter is altered
        funandgames = numpy.split(Idontsing[u], 2)#sets up so that we can alter the specific transporter end
        if h in y_cord[j]:
            if j in Reject_Type_One:
                E[j] = math.inf
                print("type one rejected")
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
                E[j] = math.inf
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
                print("transporter closed")
                print("rejected=", P[j], j)
            elif numpy.array_equal(downone, funandgames[0]):
                E[j] = 0
                print("glitch", funandgames[0], funandgames[1])
            else:
                E[j] = 0 
    elif DC1 in x_cord[j]:#this specifically will need MAJOR adjustements if we loop all instead of randomly select 
        E[j] = math.inf#energy set to infinite if particle runs into x boundaries
        print("bound", E[j], P[j])
    elif DC3 in y_cord[j]:#this isn't the best solution, as it runs through everything a couple times, however, it does have 100% accuracy
        E[j] = math.inf#energy set to infinite if particle runs into y boundaries and membrane
        print("bound", E[j], P[j])
    else:
        E[j] = 0
def func2():
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
def func3():
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
    type1low = 0
    type1high = 0
    type2low = 0
    type2high = 0
    y_cord_old = x[1].copy()
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
                break #if the particle moves into the same position as another particle, it's energy is infinite and the loop stops there
            else:
                o = 0
    if j in range(0, (k + c + g + f)) and o == 0:
        func1()
    elif j in range((k + c + g + f), (Transrand + k + c + g + f)):
        print("yay")
        u = j - (k + c + g + f)       
        funandgames = numpy.hsplit(Idontsing[u], 2)
        if funandgames[0] in upone:
            print("cool", u, j)
            funandgames[0] = downone
            if funandgames[1] in downone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
            elif funandgames[1] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
        elif funandgames[0] in downone:
            print("cool", u, j)
            funandgames[0] = upone
            if funandgames[1] in downone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
            elif funandgames[1] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == uph and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
        print("E[j]=", E[j])
        Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
        print("lol", i, j, P, Idontsing)
        for j in range(0, (k + c + g + f)):
            func1()
        for j in range((k + c + g + f), ((Transrand) + k + c + g + f)):
            func2()
        for j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):
            func3()
    elif j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):
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
            elif funandgames[0] in downone:      
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young + partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
        elif funandgames[1] in downone:
            print("cool", u, j) 
            funandgames[1] = upone
            if funandgames[0] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = epsilion + E_young
            elif funandgames[0] in upone:
                for j in range(0, (k + c + g + f)):
                    if y_cord[j] == h and x_cord[j] in Sppecial[u]:
                        E[j] = E_young - partcount
                        E_young = E[j].copy()
                E[j] = negepsilion + E_young
        print("E[j]=", E[j])
        Idontsing[u] = numpy.concatenate((funandgames[0], funandgames[1]), axis= None)
        print("lol", i, j, P, Idontsing)
        for j in range(0, (k + c + g + f)):
            func1()
        for j in range((k + c + g + f), ((Transrand) + k + c + g + f)):
            func2()
        for j in range((Transrand + k + c + g + f), ((2 * Transrand) + k + c + g + f)):
            func3()
    print("norm", i, P, E, Idontsing)
    if E[j] > Ej_old:#if the energy's lower or equal we automatically accept it
        val = numpy.random.uniform(low=0, high=1, size=1)
        if val > (math.exp((Ej_old - E[j])/(KbT))):#equation will need to be changed later(rn inaccurate)
            #T represents thermal energy, unsure how we want that represented(Kb * T)
            if j in range(0, (k + c + g + f)):
                P[j] = P_old
                E = E_old
                Idontsing = I_old
            if j in range((k + c + g + f), ((2 * Transrand) + (k + c + g + f))):
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
    type1highs[i] = type1high#remove these four lines when not doing diffussion testing
    type2highs[i] = type2high
    type1lows[i] = type1low
    type2lows[i] = type2low
    densitylow = type1low + type2low #number of particles in the lower compartment
    densityhigh = type1high + type2high #number of particles in the higher compartment
    print("type1low=", type1low/arealow, "type2low=", type2low/arealow, "lowdensity=", densitylow/arealow, file= fileout)
    print("type1high=", type1high/areahigh, "type2high=", type2high/areahigh, "highdensity=", densityhigh/areahigh, file= fileout)
fileout.close() 
fileout = open ("p.txt", "w")
for j in range(0, (k + c + g +f)):
    print(*P[j], file= fileout)
fileout.close()
pylab.title("highgraph type1") 
pylab.xlim([0, n])
pylab.ylim([0, 40])
pylab.plot(lol, type1highs)
pylab.savefig("highgrapht1.png",dpi=600) 
pylab.show() 
pylab.title("lowgraph type1") 
pylab.xlim([0, n])
pylab.ylim([0, 40])
pylab.plot(lol, type1lows)
pylab.savefig("lowgrapht1.png",dpi=600) 
pylab.show() 
pylab.title("highgraph type2") 
pylab.xlim([0, n])
pylab.ylim([0, 40])
pylab.plot(lol, type2highs)
pylab.savefig("highgrapht2.png",dpi=600) 
pylab.show()
pylab.title("lowgraph type2") 
pylab.xlim([0, n])
pylab.ylim([0, 40])
pylab.plot(lol, type2lows)
pylab.savefig("lowgrapht2.png",dpi=600) 
pylab.show()