from pyhyp import pyHyp

fileName = "merged.xyz"
fileType = "plot3d"

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": fileName,
    "fileType": fileType,
    "unattachedEdgesAreSymmetry": False,
    "outerFaceBC": "overset",
    "autoConnect": True,
    "BC": {},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 25,
    "s0": 1e-2,
    "marchDist": 60.5,
    "splay": 0.7,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1,
    "pGridRatio": -1,
    "cMax": 5,
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    "epsE": 1.0,
    "epsI": 2.0,
    "theta": 3.0,
    "volCoef": 0.25,
    "volBlend": 0.0001,
    "volSmoothIter": 100,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS("collar.cgns")
