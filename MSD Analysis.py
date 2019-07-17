import numpy
import matplotlib.pyplot as plt
import pylab
import random

n=10 #number of timesteps taken
k=3 #number of particles
d=2 #number of dimensions
#Possibly use "f_name=" to save time editing when using multiple file names
D = numpy.genfromtxt("coordinatesold.txt")
 #([row, column])
MSD_List = []
DC_List = []
fileout = open ("msdanalysis.txt", "w")
print("i  ", end = " ", file= fileout) 
for q in range(0, k):
    print("P[", q, "]", "  ", end = " ", file= fileout)#might need different spacing for all d=1,2,3
    if q == k - 1:
        print(file= fileout)
for i in range(0, 1):
    print(i, D[i], file= fileout)
for i in range(1, n):
    MSD = sum(((D[i] - D[0])**2)/k)#made up term used to calculate MSD
    MSD_List.append(MSD)
    DC = (MSD/(2 * d * i))
    DC_List.append(DC)
    print(i, D[i], end = " ", file= fileout)
    print("MSD=", MSD, "DC=", DC, file= fileout)
fileout.close()
fig, ax = plt.subplots()
ax.set(xlabel='time (s)', ylabel='Mean Squared Displacement',
       title="Mean Squared Displacement over Time($n = " + str(n) + "$ steps, $k = " + str(k) + "$ particles)")
ax.grid()
ax.plot(range(1, n), MSD_List)
plt.show()
fig.savefig("MSD.png")
fig, ax = plt.subplots()
ax.set(xlabel='time (s)', ylabel='Diffussion Constant',
       title="Diffussion Constant over Time($n = " + str(n) + "$ steps, $k = " + str(k) + "$ particles)")
ax.grid()
ax.plot(range(1, n), DC_List)
plt.show()
fig.savefig("DC.png")
