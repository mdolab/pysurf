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

outputDirectory = './'
saveRepositoryInfo(outputDirectory)

# Define input files
ransFile = '../crm_wb_L2.cgns' # specify the initial CFD grid

# Default to comm world
comm = MPI.COMM_WORLD

# Constants
in2m = 0.0254

# Converting units (data taken from DPW6 Geoemtry website)
cRef = 275.80*in2m
sRef = 594720*in2m*in2m*0.5
bRef = 2313.50*in2m
xRef = 1325.90*in2m
yRef =  468.75*in2m
zRef =  177.95*in2m

# Set guess for angle of attack
alpha0 = 2.4

# Common aerodynamic problem description and design variables
ap = AeroProblem(name='fc', alpha=alpha0,
                 mach=0.85, reynolds=5.0e6, reynoldsLength=cRef, T=310.928,
                 areaRef=sRef, chordRef=cRef, xRef=xRef, yRef=yRef, zRef=zRef,
                 evalFuncs=['cl','cd','cmy'])
CLstar = 0.5

ap.addDV('alpha')

AEROSOLVER = ADFLOW

gridFile = ransFile
CFL = 1.0
MGCYCLE = 'sg'
MGSTART = 1
useNK = False
L2Conv = 1e-7

aeroOptions = {
    # Common Parameters
    'gridFile':gridFile,
    'outputDirectory':outputDirectory,
    #'restartFile':'fc_restart_vol.cgns',
    'writeSurfaceSolution':True,
    'writeVolumeSolution':True,
    'isoSurface':{'shock':1,'vx':-0.0001},
    'solutionPrecision':'double',

    # Physics Parameters
    'equationType':'rans',
    'smoother':'dadi',
    'nsubiter':3,
    'nsubiterturb':3,
    'liftIndex':3,
    'vis2':0.250,
    #'vis4':vis4,
    #'turbulenceorder':'first order',
    #'discretization':discretization,
    #'useqcr':qcr,
    #'userotationsa':qcr,
    #'turbulenceproduction':turbulenceproduction,

    # Common Parameters
    'CFL':CFL,
    'CFLCoarse':CFL,
    'MGCycle':MGCYCLE,
    'MGStartLevel':MGSTART,
    'nCyclesCoarse':3,
    'nCycles':5, # 7500,
    'monitorvariables':['resrho','resturb', 'cl','cd','cdp','cdv','cmy', 'cpu'],
    'surfaceVariables':['vx','vy','vz','rho','P','cp','cf','cfx','cfy','cfz','blank'],
    'volumevariables':['resrho','mach','cp','resturb','blank'],

    # NK parameters
    'useNKSolver':useNK,
    'nkjacobianlag':2,
    'nkasmoverlap':2,
    'nkpcilufill':1,
    'blocksplitting':True,
    'loadbalanceiter':2,
    'nkswitchtol': 7.0e-7,
    'nkouterpreconits':5,
    'nkinnerpreconits':2,
    'nkls':'non monotone',

    # Convergence Parameters
    'L2Convergence':L2Conv,
    'L2ConvergenceCoarse':1e-4,
    'miniterationnum':250,
    'rkreset':True,

    # Overset parameters
    'nearWallDist':0.1,
    'debugzipper':True,
}

# Create solver
CFDSolver = AEROSOLVER(options=aeroOptions, comm=comm, debug=True)

'''
# Uncoment this if just want to check flooding
CFDSolver.setAeroProblem(ap)
CFDSolver.writeSolution()
'''

# Solve and evaluate
funcs = {}
CFDSolver(ap)

CFDSolver.evalFunctions(ap, funcs)

# Evaluate Gradients
#funcsSens = {}
#CFDSolver.evalFunctionsSens(ap, funcsSens)

if MPI.COMM_WORLD.rank == 0:
    print('Functions:', funcs)
    #print 'Functions Sens:', funcsSens
