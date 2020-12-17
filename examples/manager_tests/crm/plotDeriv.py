# IMPORTS
import numpy as np
import pickle
import matplotlib.pyplot as plt

# Load file with results
with open("results.pickle", "r") as fid:
    results = pickle.load(fid)

# Split data
Z = results[0, :]
Y = results[1, :]
dYdZ = results[2, :]

# Compute arrow coordinates
dZ = Z[1] - Z[0]
Xarrow = np.ones(len(Z)) * dZ
Yarrow = dYdZ[:] * dZ

from matplotlib import patches as mpatches


def arrow(self, x, y, dx, dy, **kwargs):
    kwargs.setdefault("arrowstyle", "simple, head_width=10, head_length=10")
    kwargs.setdefault("fc", "black")
    x = self.convert_xunits(x)
    y = self.convert_yunits(y)
    dx = self.convert_xunits(dx)
    dy = self.convert_yunits(dy)
    posA = x, y
    posB = x + dx, y + dy
    a = mpatches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
    self.add_artist(a)
    return a


# Plot
fig = plt.figure()
plt.plot(Y, Z, "o-")
ax = plt.axes()
for ii in range(len(Z)):
    arrow(ax, Y[ii], Z[ii], Yarrow[ii], Xarrow[ii], fc="r", ec="r")
plt.axis("equal")
plt.xlabel("Y")
plt.ylabel("Z")

plt.show()

fig.savefig("Z_sweep.png", dpi=300)
