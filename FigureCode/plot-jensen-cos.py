from __future__ import division  # 2.2+-only
import matplotlib.pyplot as plt     # For Matlab-like plotting
import numpy as np
import scipy as sp
import math                         # For tan() function

# Initial conditions of Benchmark Case 1
Uinf = 8.10     # (m/s), for the unfirom wind field case given in (Jensen, 1983) pg6.
r = 20          # (m), same as above, radius of rotor
Alpha = 0.1     # Entrainment Constant, (Jensen, 1983) pg6.
x = [40, 100]   # Two distances enumerated in (Jensen, 1983) pg6.
indexX = 1      # Which x to use
Ratio = [16, 10, 6]     # Three ratios appearing in (Jensen, 1983) pg7.
VelRotorFlat = [0 for i in range(3)]  # initialize with zeros
VelRtrCos = [[0 for j in range(20)]
             for i in range(3)]   # initialize with zeros
VStar = [0,0,0]  #initialize with zeros
VStarCos = [[0 for j in range(20)] for i in range(3)]  # initialize with zeros

Theta = math.degrees(math.atan((50/3)/(500/3)))

for i in range(0, 3):
    #print "Ratio[",i,"] = ", Ratio[i] # debug print statement
    rKnot = x[indexX] / Ratio[i] # from x/ro = 16, top graph, (Jensen, 1983) pg7.
    VelRotorFlat[i] = Uinf *  (1 - ((2/3) * ((rKnot / (rKnot + (Alpha * x[indexX] ))) ** 2)))  # Taken from pg6 of (Jensen, 1983).
    VStar[i] = (VelRotorFlat[i] / Uinf)
    # Now do our cosine curve
    for Angle in range(0, 20):
        # Taken from pg8 of (Jensen, 1983).
        f = (1 + math.cos(math.radians(9 * Angle)))/2
        VelRtrCos[i][Angle] = Uinf * (1 - ((2/3) * (((f * rKnot) / (rKnot + (Alpha * x[indexX]))) ** 2)))
        VStarCos[i][Angle] = (VelRtrCos[i][Angle] / Uinf)
        #print VStarCos[i][Angle]
    


#print "x = ", x[indexX]
#print "Ratio = ", Ratio[indexR]
#print "V in the wake is: ", VelRotor[i]
#print "v/u is", (VelRotor[i] / Uinf)

# For checking what's in our Cos array
#for i in range(len(VelRtrCos)):
# for j in range(len(VelRtrCos[1])):
#  print VelRtrCos[1][j]

xRange = [-30, -20, -10, -(Theta+.0001), -Theta, 0, Theta,(Theta+.0001), 10, 20, 30]     # Taken from (Jensen, 1983) pg7.
xRangeCos = np.linspace(-30, 30, num=61)
yValsFlat1 = [1, 1, 1, 1, VStar[0], VStar[0], VStar[0], 1, 1, 1, 1]
yValsFlat2 = [1, 1, 1, 1, VStar[1], VStar[1], VStar[1], 1, 1, 1, 1]
yValsFlat3 = [1, 1, 1, 1, VStar[2], VStar[2], VStar[2], 1, 1, 1, 1]
yValsFlat = [yValsFlat1, yValsFlat2, yValsFlat3]
xRanLast = len(xRange) - 1
yRange1 = [0.8, 0.9, 1.0, 1.1]
yRange2 = [0.8, 0.9, 1.0]
yRange3 = [0.7, 0.8, 0.9, 1.0]
yRange = [yRange1, yRange2, yRange3]
yRanLast = len(xRange) - 1

plt.close('all')

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(xRange, yValsFlat[0])
# Flat Distribution
ax1.axis([xRange[0], xRange[xRanLast], 0.8, 1.1])
# Cos Distribution
ax1.plot(xRangeCos, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     VStarCos[0][19], VStarCos[0][18], VStarCos[0][17], VStarCos[0][16], VStarCos[0][15],
                     VStarCos[0][14], VStarCos[0][13], VStarCos[0][12], VStarCos[0][11], VStarCos[0][10],
                     VStarCos[0][9], VStarCos[0][8], VStarCos[0][7], VStarCos[0][6], VStarCos[0][5],
                     VStarCos[0][4], VStarCos[0][3], VStarCos[0][2], VStarCos[0][1], VStarCos[0][0],
                     VStarCos[0][1], VStarCos[0][2], VStarCos[0][3], VStarCos[0][4], VStarCos[0][5],
                     VStarCos[0][6], VStarCos[0][7], VStarCos[0][8], VStarCos[0][9], VStarCos[0][10],
                     VStarCos[0][11], VStarCos[0][12], VStarCos[0][13], VStarCos[0][14], VStarCos[0][15],
                     VStarCos[0][16], VStarCos[0][17], VStarCos[0][18], VStarCos[0][19], 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
plt.yticks(np.arange(min(yRange[0]), max(yRange[0]), 0.1))

# Flat Distribution
ax2.plot(xRange, yValsFlat[1])
# Cos Distribution
ax2.plot(xRangeCos, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,
                     VStarCos[1][19], VStarCos[1][18], VStarCos[1][17], VStarCos[1][16], VStarCos[1][15],
                     VStarCos[1][14], VStarCos[1][13], VStarCos[1][12], VStarCos[1][11], VStarCos[1][10],
                     VStarCos[1][9], VStarCos[1][8], VStarCos[1][7], VStarCos[1][6], VStarCos[1][5],
                     VStarCos[1][4], VStarCos[1][3], VStarCos[1][2], VStarCos[1][1], VStarCos[1][0],
                     VStarCos[1][1], VStarCos[1][2], VStarCos[1][3], VStarCos[1][4], VStarCos[1][5],
                     VStarCos[1][6], VStarCos[1][7], VStarCos[1][8], VStarCos[1][9], VStarCos[1][10],
                     VStarCos[1][11], VStarCos[1][12], VStarCos[1][13], VStarCos[1][14], VStarCos[1][15],
                     VStarCos[1][16], VStarCos[1][17], VStarCos[1][18], VStarCos[1][19],1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
ax2.axis([xRange[0], xRange[xRanLast], 0.8, 1.1])
plt.yticks(np.arange(min(yRange[1]), max(yRange[1]), 0.1))

# Flat Distribution
ax3.plot(xRange, yValsFlat[2])
# Cos Distribution
ax3.plot(xRangeCos, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     VStarCos[2][19], VStarCos[2][18], VStarCos[2][17], VStarCos[2][16], VStarCos[2][15],
                     VStarCos[2][14], VStarCos[2][13], VStarCos[2][12], VStarCos[2][11], VStarCos[2][10],
                     VStarCos[2][9], VStarCos[2][8], VStarCos[2][7], VStarCos[2][6], VStarCos[2][5],
                     VStarCos[2][4], VStarCos[2][3], VStarCos[2][2], VStarCos[2][1], VStarCos[2][0],
                     VStarCos[2][1], VStarCos[2][2], VStarCos[2][3], VStarCos[2][4], VStarCos[2][5],
                     VStarCos[2][6], VStarCos[2][7], VStarCos[2][8], VStarCos[2][9], VStarCos[2][10],
                     VStarCos[2][11], VStarCos[2][12], VStarCos[2][13], VStarCos[2][14], VStarCos[2][15],
                     VStarCos[2][16], VStarCos[2][17], VStarCos[2][18], VStarCos[2][19], 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
ax3.axis([xRange[0], xRange[xRanLast], 0.6, 1.1])

# To get the y and tick marks like (Jensen 83)
plt.yticks(np.arange(min(yRange[2]), max(yRange[2]), 0.1))

plt.ylabel('v/u')
plt.xlabel(r'$\Delta$ ' + r'$\theta$ ' + '[deg.]')

# Adjust the subplot layout, because the names might take more space than usual"
f.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)

plt.show()
