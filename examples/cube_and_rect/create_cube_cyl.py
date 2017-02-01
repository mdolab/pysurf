"""
This is an example script to show the entire process of creating overset
meshes using the MDOlab's tools. Here, we have a cube that is intersecting
a rectangular prism with one of the cube's edges being a symmetry plane.

|
|---------------+
|               |
|               +---------------------+
|               |                     |
|               +---------------------+
|               |          rect
|---------------+
|         cube
sym

We defined functions to call individual tools from the group to create
these meshes. Each will be explained in its docstring.

John Jasa 2017-01

"""

# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import os
from pyhyp import pyHyp

def extrude_cube_volume_mesh():
    """
    First we need to create a primary volume mesh for the cube.
    We already have a structured surface mesh in `cube.cgns`.
    Next, we use pyHyp to hyperbollicaly extrude a volume mesh.
    Note that we need to set the symmetry boundary condition in the options.
    """

    # Input filename for pyHyp
    fileName = 'cube.cgns'

    options= {
        # ---------------------------
        #        Input Parameters
        # ---------------------------
        'inputFile':fileName,
        'fileType':'CGNS',
        'unattachedEdgesAreSymmetry':True,
        'outerFaceBC':'overset',
        'autoConnect':True,
        'families':'wall',

        # ---------------------------
        #        Grid Parameters
        # ---------------------------
        'N': 9,
        's0': 1.e-1,
        'marchDist':1.2,
        'splay':0.2,

        # ---------------------------
        #   Pseudo Grid Parameters
        # ---------------------------
        'ps0': -1,
        'pGridRatio': -1,
        'cMax': 50,

        # ---------------------------
        #   Smoothing parameters
        # ---------------------------
        'epsE': 0.10,
        'epsI': 2.0,
        'theta': 3.0,
        'volCoef': 0.25,
        'volBlend': 0.01,
        'volSmoothIter': 400,
        }

    hyp = pyHyp(options=options)
    hyp.run()
    hyp.writeCGNS('cube_vol.cgns')


def extrude_rect_volume_mesh():
    """
    Next we need to create a primary volume mesh for the rectangular prism.
    We already have a structured surface mesh in `rect.cgns`.
    Next, we use pyHyp to hyperbollicaly extrude a volume mesh.
    Note that we need to set the symmetry boundary condition in the options.
    """
    fileName = 'rect.cgns'

    options= {
        # ---------------------------
        #        Input Parameters
        # ---------------------------
        'inputFile':fileName,
        'fileType':'CGNS',
        'unattachedEdgesAreSymmetry':True,
        'outerFaceBC':'overset',
        'autoConnect':True,
        'BC':{},
        'families':'wall',

        # ---------------------------
        #        Grid Parameters
        # ---------------------------
        'N': 9,
        's0': 1.e-1,
        'marchDist':1.2,
        'splay':0.2,

        # ---------------------------
        #   Pseudo Grid Parameters
        # ---------------------------
        'ps0': -1,
        'pGridRatio': -1,
        'cMax': 50,

        # ---------------------------
        #   Smoothing parameters
        # ---------------------------
        'epsE': .10,
        'epsI': 1.0,
        'theta': 3.0,
        'volCoef': 0.25,
        'volBlend': 0.01,
        'volSmoothIter': 400,
        }

    hyp = pyHyp(options=options)
    hyp.run()
    hyp.writeCGNS('rect_vol.cgns')


def march_surface_meshes():
    """
    Now we will march the surface meshes on both the cube and the rectangle
    using the intersection curve as the starting curve.
    We do this using hypsurf, the included surface meshing tool.
    """

    # Set communicator
    comm = MPI.COMM_WORLD

    # Load the unstructured cube and rectangular prism cgns files
    cube = pysurf.TSurfGeometry('cube_uns.cgns', comm)
    rect = pysurf.TSurfGeometry('rect_uns.cgns', comm)

    # Intersect the cube and rect and obtain the first (and only) intersection curve
    curve = cube.intersect(rect)[0]

    # Split curves based on sharpness
    splitcurves = pysurf.tsurf_tools.split_curve_single(curve, 'intersections', criteria='sharpness')

    # Remesh each individual edge of the intersection curve
    for curveName in splitcurves:
        splitcurves[curveName] = splitcurves[curveName].remesh(nNewNodes=21)

    # Merge the four curves back together into a single curve
    curve = pysurf.tsurf_tools.merge_curves(splitcurves,'intersection')

    # Export intersection curve
    curve.export_tecplot('curve')

    # Create curve list based on imported curves
    # !!! Make sure to call `extract_curves.py` before running this script
    curves = []
    curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_000.plt'))
    curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_001.plt'))
    curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_002.plt'))
    curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_003.plt'))

    # Create an empty list in which we'll store the long (longitudinal) edges of
    # the rect
    long_curves = []
    counter = 0

    # Examine each of the imported curves
    for ext_curve in curves:

        # Split these curves based on sharpness to get all edges of the rect
        split_curve = pysurf.tsurf_tools.split_curve_single(ext_curve[ext_curve.keys()[0]], 'int', criteria='sharpness')

        # Loop over these split curves
        for name in split_curve:

            # Give the curves new names so they do not conflict with each other
            split_curve[name].name = 'int_'+'{}'.format(counter).zfill(3)
            counter += 1

            # Extract the curve points and compute the length
            pts = split_curve[name].extract_points()
            length = pts[0, 2] - pts[-1, 2]

            # Hardcoded logic here based on the length of edges
            if np.abs(length) > 1:

                # Flip the long curve if it's facing the 'wrong' direction
                if length > 0:
                    split_curve[name].flip()

                # Add the long curve to the list
                long_curves.append(split_curve[name])

    # Create a list of guideCurves based on the extracted long curves
    # Note here we need the strings, not the curve objects
    guideCurves = []
    for ext_curve in long_curves:
        rect.add_curve(ext_curve)
        guideCurves.append(ext_curve.name)

    # Add intersection curve to each component
    rect.add_curve(curve)
    cube.add_curve(curve)

    # Function to actually call hypsurf given inputs for the cube and rect
    def generateMesh(curve, geom, options, output_name):

        # Set guideCurves for the rect case
        if output_name == 'rect.xyz':
            options['remesh'] = True
            options['guideCurves'] = guideCurves
        else:
            options['remesh'] = False
            options['guideCurves'] = []

        mesh = pysurf.hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

        mesh.createMesh()

        mesh.exportPlot3d(output_name)

    # GENERATE RECT MESH
    print 'Generating rect mesh'

    # Set options
    options = {
        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.03,
        'numLayers' : 17,
        'extension' : 3.5,
        'epsE0' : 4.5,
        'theta' : -0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'nuArea' : 0.16,
        'numAreaPasses' : 20,
        'sigmaSplay' : 0.3,
        'cMax' : 10000.0,
        'ratioGuess' : 1.5,
    }

    # Flip the curve so the normal orientation is correct.
    # This is very important so the marching function knows which direction
    # to march on the surface.
    curve.flip()

    # Call meshing function
    generateMesh(curve, rect, options, 'rect.xyz')

    # GENERATE CUBE MESH
    print 'Generating cube mesh'

    # Options
    options['dStart'] = 0.02
    options['numLayers'] = 17
    options['extension'] = 2.5

    # Flip the curve back to its original orientation.
    curve.flip()

    # Call meshing function
    generateMesh(curve, cube, options, 'cube.xyz')


def merge_surface_meshes():
    """
    This function takes the two blocks that define the collar mesh
    and joins them so that we can run pyHyp
    """

    pysurf.plot3d_interface.merge_plot3d(['cube.xyz', 'rect.xyz'], [[1,0,0], [0,1,0]])

    mergedGrid = pysurf.plot3d_interface.read_plot3d('merged.xyz',3)

    mergedGrid.remove_curves()

    pysurf.plot3d_interface.export_plot3d(mergedGrid, 'merged.xyz')


def run_pyhyp_for_collar():
    """
    Now that we have the merged surface mesh for the collar between the
    cube and rect, we need to run pyHyp to extrude a volume mesh for the collar.
    """

    print
    print 'Now running pyHyp on the merged mesh'
    print

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
        'N': 35,
        's0': 1e-2,
        'marchDist':1.3,
        'splay':0.7,
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


def create_OCart_mesh_and_merge():
    """
    Now that we have two primary meshes and the collar mesh, we only need a
    background mesh to complete our overset mesh system.
    To create this background mesh, we use cgns_utils to create an OCart mesh,
    which has a uniformly sized cartesian mesh overlapping the existing bodies,
    then calls pyHyp to extrude an O mesh around this block.
    """

    print
    print "Now creating OCart mesh"
    print

    import subprocess

    subprocess.call(["cgns_utils combine cube_vol.cgns rect_vol.cgns collar.cgns half_full.cgns"], shell=True)

    # Now create the background mesh
    # positional arguments:
    #  gridFile    Name of input CGNS file
    #  dh          Uniform cartesian spacing size
    #  hExtra      Extension in "O" dimension
    #  nExtra      Number of nodes to use for extension
    #  sym         Normal for possible sym plane
    #  mgcycle     Minimum MG cycle to enforce
    #  outFile     Name of output CGNS file
    subprocess.call(["cgns_utils simpleOCart half_full.cgns 0.1 5. 9 z 1 background.cgns"], shell=True)

    # Now combine everything in a single file
    subprocess.call(["cgns_utils combine half_full.cgns background.cgns full.cgns"], shell=True)

    # Use cgns_utils to merge contiguous blocks
    subprocess.call(["cgns_utils merge full.cgns"], shell=True)

    # Check block-to-block connectivities
    subprocess.call(["cgns_utils connect full.cgns"], shell=True)


def run_ADflow_to_check_connections():
    """
    Lastly, we run ADflow to check the domain connectivity between the overset
    meshes. Note that it doesn't make sense to actually run CFD on these
    bodies, but this serves as a simple example case to create a collar mesh.
    """

    # ======================================================================
    #         Import modules
    # ======================================================================
    import numpy
    import argparse
    from mpi4py import MPI
    from baseclasses import AeroProblem
    from tripan import TRIPAN
    from adflow import ADFLOW
    from repostate import saveRepositoryInfo

    # ======================================================================
    #         Input Information
    # ======================================================================

    outputDirectory = './'
    saveRepositoryInfo(outputDirectory)

    # Default to comm world
    comm = MPI.COMM_WORLD

    # Common aerodynamic problem description and design variables
    ap = AeroProblem(name='fc', alpha=1.4, mach=0.1, altitude=10000,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cd'])

    AEROSOLVER = ADFLOW

    CFL = 1.0
    MGCYCLE = 'sg'
    MGSTART = 1
    useNK = False

    aeroOptions = {
        # Common Parameters
        'gridFile':'full.cgns',
        'outputDirectory':outputDirectory,

        # Physics Parameters
        'equationType':'rans',
        'smoother':'dadi',
        'liftIndex':2,

        # Common Parameters
        'CFL':CFL,
        'CFLCoarse':CFL,
        'MGCycle':MGCYCLE,
        'MGStartLevel':MGSTART,
        'nCyclesCoarse':250,
        'nCycles':1000,
        'monitorvariables':['resrho','cl','cd'],
        'volumevariables':['resrho','blank'],
        'surfacevariables':['cp','vx', 'vy','vz', 'mach','blank'],
        'nearWallDist':0.1,
        'nsubiterturb':3,
        'useNKSolver':useNK,

        # Debugging parameters
        'debugzipper':True,
    }

    # Create solver
    CFDSolver = AEROSOLVER(options=aeroOptions, comm=comm, debug=True)

    # Uncoment this if just want to check flooding
    CFDSolver.setAeroProblem(ap)
    CFDSolver.writeSolution()

# Here we actually run all the functions we just defined.
extrude_cube_volume_mesh()
extrude_rect_volume_mesh()
march_surface_meshes()
merge_surface_meshes()
run_pyhyp_for_collar()
create_OCart_mesh_and_merge()
run_ADflow_to_check_connections()
