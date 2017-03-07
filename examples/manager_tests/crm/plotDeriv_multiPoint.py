# IMPORTS
from __future__ import division
import numpy as np
import pickle
import matplotlib.pyplot as plt

# Give nodal positions
numNodes = 11

# WING POSITIONS
deltaZ = np.linspace(0.000001, 3.0, 11)

def load_data(nodeID):

    # Load file with results
    with open('results_%03d.pickle'%(nodeID),'r') as fid:
        results = pickle.load(fid)

    # Split data
    Y = results[1][0,:]
    dYdZ = results[1][1,:]
    Z = np.ones(len(Y))*results[0]

    # Compute arrow coordinates
    dZ = deltaZ[1]-deltaZ[0]
    Xarrow = np.ones(len(Y))*dZ
    Yarrow = dYdZ[:]*dZ

    return Y, Z, Xarrow, Yarrow

# LOAD DATA FOR ALL NODAL POSITIONS

# Initialize lists to store data
Y_tot = []
Z_tot = []
Xarrow_tot = []
Yarrow_tot = []

# Loop for every nodal position
for nodeID in range(numNodes):

    # Load data from pickle file
    Y, Z, Xarrow, Yarrow = load_data(nodeID)

    # Store information to the lists
    Y_tot.append(Y)
    Z_tot.append(Z)
    Xarrow_tot.append(Xarrow)
    Yarrow_tot.append(Yarrow)

# DEFINE FUNCTION TO PLOT ARROWS

from matplotlib import patches as mpatches

def arrow(self, x, y, dx, dy, **kwargs):
    kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
    kwargs.setdefault('fc', 'black')
    x = self.convert_xunits(x)
    y = self.convert_yunits(y)
    dx = self.convert_xunits(dx)
    dy = self.convert_yunits(dy)
    posA = x, y
    posB = x+dx, y+dy
    a = mpatches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
    self.add_artist(a)
    return a

# PLOT
fig = plt.figure()

for nodeID in range(numNodes):

    Y = Y_tot[nodeID]
    Z = Z_tot[nodeID]
    Xarrow = Xarrow_tot[nodeID]
    Yarrow = Yarrow_tot[nodeID]

    plt.plot(Y,Z,'o',label='j=%02d'%(nodeID))
    ax = plt.axes()
    for ii in range(len(Z)):
        arrow(ax, Y[ii], Z[ii], Yarrow[ii], Xarrow[ii], fc='r', ec='r')


plt.axis('equal')
plt.xlabel('Y')
plt.ylabel('Z')
#plt.legend(loc='best')

plt.show()

fig.savefig('Z_sweep.png',dpi=300)
