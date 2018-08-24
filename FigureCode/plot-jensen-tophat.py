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
VelRotor = [0,0,0]  #initialize with zeros
VStar = [0,0,0]  #initialize with zeros

Theta = math.degrees(math.atan((50/3)/(500/3)))

for i in range(0, 3):
    #print "Ratio[",i,"] = ", Ratio[i] # debug print statement
    rKnot = x[indexX] / Ratio[i] # from x/ro = 16, top graph, (Jensen, 1983) pg7.
    VelRotor[i] = Uinf *  (1 - ((2/3) * ((rKnot / (rKnot + (Alpha * x[indexX] ))) ** 2)))  # Taken from pg6 of (Jensen, 1983).
    VStar[i] = (VelRotor[i] / Uinf)

#print "x = ", x[indexX]
#print "Ratio = ", Ratio[indexR]
#print "V in the wake is: ", VelRotor[i]
#print "v/u is", (VelRotor[i] / Uinf)

xRange = [-30, -20, -10, -(Theta+.0001), -Theta, 0, Theta,(Theta+.0001), 10, 20, 30]     # Taken from (Jensen, 1983) pg7.
xRanLast = len(xRange) - 1
yRange1 = [0.8, 0.9, 1.0, 1.1]
yRange2 = [0.8, 0.9, 1.0]
yRange3 = [0.7, 0.8, 0.9, 1.0]
yRange = [yRange1, yRange2, yRange3]
yRanLast = len(xRange) - 1

plt.close('all')

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(xRange, [1,1,1,1, VStar[0],VStar[0],VStar[0],1, 1,1,1])
ax1.axis([xRange[0], xRange[xRanLast], 0.8, 1.1])
ax2.plot(xRange, [1,1,1,1, VStar[1],VStar[1],VStar[1],1, 1,1,1])
ax2.axis([xRange[0], xRange[xRanLast], 0.8, 1.1])
ax3.plot(xRange, [1,1,1,1, VStar[2],VStar[2],VStar[2],1, 1,1,1])
ax3.axis([xRange[0], xRange[xRanLast], 0.6, 1.1])
plt.ylabel('v/u')
plt.xlabel('Theta(deg)')

# Adjust the subplot layout, because the names might take more space than usual"
f.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)

plt.show()
