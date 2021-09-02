# IMPORTS
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches

# INITIALIZE FIGURE
fig = plt.figure()

# DEFINE GENERAL PLOTTING FUNCTION
def plotLevel(level, maxLevel):

    # Load file with results
    with open("results_L%d.pickle" % level, "r") as fid:
        results = pickle.load(fid)

    # Split data
    Z = results[0, :]
    Y = results[1, :]
    dYdZ = results[2, :]

    # Compute arrow coordinates
    dZ = Z[1] - Z[0]
    Xarrow = np.ones(len(Z)) * dZ
    Yarrow = dYdZ[:] * dZ

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

    # Compute transparency factor
    alpha = 0.2 + 0.8 * (1 - level / maxLevel)

    # Plot
    plt.plot(Y, Z, "ob", alpha=alpha, markeredgewidth=0.0)
    ax = plt.axes()
    for ii in range(len(Z)):
        arrow(ax, Y[ii], Z[ii], Yarrow[ii], Xarrow[ii], fc="r", ec="r", alpha=alpha)
        plt.axis("equal")
        plt.xlabel("Y")
        plt.ylabel("Z")


# RUN FUNCTION FOR EVERY LEVEL
maxLevel = 3
levels = range(maxLevel + 1)

for level in levels:
    plotLevel(maxLevel - level, maxLevel)

plt.show()

fig.savefig("cylinder.png", dpi=300)
