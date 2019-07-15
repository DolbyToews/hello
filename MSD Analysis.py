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
for i in range(0, 1):
    print(i, D[i], file= fileout)
#print("MSD=", numpy.mean(MSD), file=fileout)
#print("diffussion constant=", numpy.mean(MSD)/(2 * d * i), file=fileout)
#print(file= fileout)#this forces it to make a new line after the row is complete
for i in range(1, n):
    #for x in range(1, (2 * k) + 1):
        #MSD = ((D[i, x] - D[0,x]))
    MSD = sum(((D[i] - D[0])**2)/k)#made up term used to calculate MSD
    MSD_List.append(MSD)
    DC = (MSD/(2 * d * i))
    DC_List.append(DC)
    print(i, D[i], end = " ", file= fileout)
    print("MSD=", MSD, "DC=", DC, file= fileout)
fileout.close()
fig, ax = plt.subplots()
ax.set(xlabel='time (s)', ylabel='MSD',
       title='Means Squared Displacement over Time')
ax.grid()
ax.plot(range(1, n), MSD_List)
plt.show()
fig.savefig("MSD.png")
fig, ax = plt.subplots()
ax.set(xlabel='time (s)', ylabel='DC',
       title='Diffussion Constant modelled over Time')
ax.grid()
ax.plot(range(1, n), DC_List)
plt.show()
fig.savefig("DC.png")