# IMPORTS
import numpy as np
import pickle
import matplotlib.pyplot as plt

# Give nodal positions
i_node = 190
j_list = [0, 10, 20, 30, 40, 48]

def load_data(i_node,j_node):

    # Load file with results
    with open('results_%03d_%03d.pickle'%(i_node,j_node),'r') as fid:
        results = pickle.load(fid)

    # Split data
    Z = results[0,:]
    Y = results[1,:]
    dYdZ = results[2,:]

    # Compute arrow coordinates
    dZ = Z[1]-Z[0]
    Xarrow = np.ones(len(Z))*dZ
    Yarrow = dYdZ[:]*dZ

    return Y, Z, Xarrow, Yarrow

# LOAD DATA FOR ALL NODAL POSITIONS

# Initialize lists to store data
Y_tot = []
Z_tot = []
Xarrow_tot = []
Yarrow_tot = []

# Loop for every nodal position
for j_node in j_list:

    # Load data from pickle file
    Y, Z, Xarrow, Yarrow = load_data(i_node,j_node)

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

for jj in range(len(j_list)):

    Y = Y_tot[jj]
    Z = Z_tot[jj]
    Xarrow = Xarrow_tot[jj]
    Yarrow = Yarrow_tot[jj]

    plt.plot(Y,Z,'o-',label='j=%02d'%(j_list[jj]+1))
    ax = plt.axes()
    for ii in range(len(Z)):
        arrow(ax, Y[ii], Z[ii], Yarrow[ii], Xarrow[ii], fc='r', ec='r')


plt.axis('equal')
plt.xlabel('Y')
plt.ylabel('Z')
plt.legend(loc='best')

plt.show()

fig.savefig('Z_sweep.png',dpi=300)
