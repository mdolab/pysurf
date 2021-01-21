# ======================================================================
#         Import modules
# ======================================================================
import numpy
import argparse
from mpi4py import MPI
from baseclasses import *
from tripan import TRIPAN
from adflow import ADFLOW
from pywarp import *
from pygeo import *
from pyspline import *
from repostate import *

# ======================================================================
#         Input Information
# ======================================================================

outputDirectory = "./"
saveRepositoryInfo(outputDirectory)

# Define input files
ransFile = "../crm_wb_L2.cgns"  # specify the initial CFD grid

# Default to comm world
comm = MPI.COMM_WORLD

# Common aerodynamic problem description and design variables
ap = AeroProblem(name="fc", alpha=1.4, mach=0.80, altitude=10000, areaRef=45.5, chordRef=3.25, evalFuncs=["cl", "cd"])
ap.addDV("alpha")


AEROSOLVER = ADFLOW

gridFile = ransFile
CFL = 1.0
MGCYCLE = "sg"
MGSTART = 1
useNK = False

aeroOptions = {
    # Common Parameters
    "gridFile": gridFile,
    "outputDirectory": outputDirectory,
    #'restartFile':'fc_restart_vol.cgns',
    # Physics Parameters
    "equationType": "rans",
    "smoother": "dadi",
    "liftIndex": 2,
    # Common Parameters
    "CFL": CFL,
    "CFLCoarse": CFL,
    "MGCycle": MGCYCLE,
    "MGStartLevel": MGSTART,
    "nCyclesCoarse": 250,
    "nCycles": 1000,
    "monitorvariables": ["resrho", "cl", "cd"],
    "volumevariables": ["resrho", "blank"],
    "surfacevariables": ["cp", "vx", "vy", "vz", "mach", "blank"],
    "nearWallDist": 0.05,
    "nsubiterturb": 3,
    "useNKSolver": useNK,
    # Convergence Parameters
    "L2Convergence": 1e-6,
    "L2ConvergenceCoarse": 1e-2,
    #'miniterationnum':100,
    # Adjoint Parameters
    "adjointL2Convergence": 1e-8,
    "ADPC": False,
    "adjointMaxIter": 500,
    "adjointSubspaceSize": 150,
    "ILUFill": 2,
    "ASMOverlap": 1,
    "outerPreconIts": 3,
    # Debugging parameters
    "debugzipper": False,
    "writeTecplotSurfaceSolution": True,
    #'overlapfactor':0.99,
}

# Create solver
CFDSolver = AEROSOLVER(options=aeroOptions, comm=comm, debug=False)

# Uncomment this if just want to check flooding
CFDSolver.setAeroProblem(ap)
CFDSolver.writeSolution()

"""
# Solve and evaluate
funcs = {}
CFDSolver(ap)

CFDSolver.evalFunctions(ap, funcs)

# Evaluate Gradients
#funcsSens = {}
#CFDSolver.evalFunctionsSens(ap, funcsSens)

if MPI.COMM_WORLD.rank == 0:
    print 'Functions:', funcs
    #print 'Functions Sens:', funcsSens
"""
