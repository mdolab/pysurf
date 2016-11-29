from pyhyp import pyHyp

fileName = 'merged.xyz'
fileType = 'plot3d'

options= {

    # ---------------------------
    #        Input Parameters
    # ---------------------------
    'inputFile':fileName,
    'fileType':fileType,
    'unattachedEdgesAreSymmetry':False,
    'outerFaceBC':'overset',
    'autoConnect':True,
    'BC':{},
    'families':'wall',

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 97, 
    's0': 1e-6,
    'marchDist':4.5,
    'splay':0.2,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0': -1,
    'pGridRatio': -1,
    'cMax': 5,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 5.0,
    'epsI': 10.0,
    'theta': 3.0,
    'volCoef': 0.25,
    'volBlend': 0.01,
    'volSmoothIter': 400,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('collar.cgns')
