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
from pywarpustruct import USMesh
from scipy.spatial import cKDTree
import subprocess

def extrude_cube_volume_mesh():
    """
    First we need to create a primary volume mesh for the cube.
    We already have a structured surface mesh in `cube_struct.cgns`.
    Next, we use pyHyp to hyperbollicaly extrude a volume mesh.
    Note that we need to set the symmetry boundary condition in the options.
    """

    # Input filename for pyHyp
    fileName = '../inputs/cube_struct.cgns'

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
    We already have a structured surface mesh in `rect_struct.cgns`.
    Next, we use pyHyp to hyperbollicaly extrude a volume mesh.
    Note that we need to set the symmetry boundary condition in the options.
    """
    fileName = '../inputs/rect_struct.cgns'

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

    print rectTranslation
    # Translate wing primary mesh
    subprocess.call(["cp rect_vol.cgns rect_vol_temp.cgns"], shell=True)
    subprocess.call(["cgns_utils translate rect_vol_temp.cgns {} {} {}".format(rectTranslation[0], rectTranslation[1], rectTranslation[2])], shell=True)


def march_surface_meshes():
    """
    Now we will march the surface meshes on both the cube and the rectangle
    using the intersection curve as the starting curve.
    We do this using hypsurf, the included surface meshing tool.
    """

    # Load components
    comp1 = pysurf.TSurfGeometry('../inputs/cube_uns.cgns',['geom'])
    comp2 = pysurf.TSurfGeometry('../inputs/rect_uns.cgns',['geom'])

    name1 = comp1.name
    name2 = comp2.name

    #name1 = 'wing'
    #name2 = 'body'

    #comp1.rename(name1)
    #comp2.rename(name2)

    # ADDING GUIDE CURVES
    # Create curve dictionary based on imported curves
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
            pts = split_curve[name].get_points()
            length = pts[2, 0] - pts[2, -1]

            # Hardcoded logic here based on the length of edges
            if np.abs(length) > 1:

                # Flip the long curve if it's facing the 'wrong' direction
                if length > 0:
                    split_curve[name].flip()

                # Add the long curve to the list
                long_curves.append(split_curve[name])

    # Create a list of guideCurves based on the extracted long curves
    # Note here we need the strings, not the curve objects.
    # We use the same loop to add the guide curves to the rectangle object
    guideCurves = []
    for ext_curve in long_curves:
        print ext_curve.name
        comp2.add_curve(ext_curve)
        if ext_curve.name != 'int_011':
            guideCurves.append(ext_curve.name)

    # Rotate the rectangle in 5 degrees
    comp2.rotate(5,2)

    # Create manager object and add the geometry objects to it
    manager0 = pysurf.Manager()
    manager0.add_geometry(comp1)
    manager0.add_geometry(comp2)

    distTol = 1e-7

    # Set up integer to export different meshes
    mesh_pass = 0

    #======================================================
    # FORWARD PASS

    def forward_pass(manager):

        '''
        This function will apply all geometry operations to the given manager.
        '''

        # INTERSECT

        # Call intersection function
        intCurveNames = manager.intersect(distTol=distTol)
        intCurveName = intCurveNames[0]

        manager.intCurves[intCurveName].shift_end_nodes(criteria='maxX')

        # REMESH
        optionsDict = {
            'nNewNodes':81,
            'spacing':'linear',
            'initialSpacing':0.005,
            'finalSpacing':0.005,
        }
        remeshedCurveName = manager.remesh_intCurve(intCurveName,optionsDict)

        manager.intCurves[remeshedCurveName].export_tecplot(remeshedCurveName)

        # MARCH SURFACE MESHES
        meshName = 'mesh'

        options_rect = {

            'bc1' : 'curve:int_011',
            'bc2' : 'curve:int_011',
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
            'guideCurves':guideCurves,

        }

        options_cube = {

            'bc1' : 'continuous',
            'bc2' : 'continuous',
            'dStart' : 0.02,
            'numLayers' : 17,
            'extension' : 2.5,
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

        meshName = 'mesh'
        meshNames = manager.march_intCurve_surfaceMesh(remeshedCurveName, options0=options_cube, options1=options_rect, meshName=meshName)

        # EXPORT
        for meshName in meshNames:
            print meshName+'_'+str(mesh_pass)+'.xyz'
            manager.meshes[meshName].exportPlot3d(meshName+'_'+str(mesh_pass)+'.xyz')

        return meshNames

    # END OF forward_pass
    #======================================================

    # Call the forward pass function to the original manager
    meshNames = forward_pass(manager0)
    mesh_pass = mesh_pass + 1


def merge_surface_meshes():
    """
    This function takes the two blocks that define the collar mesh
    and joins them so that we can run pyHyp
    """

    pysurf.plot3d_interface.merge_plot3d(['mesh_0_0.xyz', 'mesh_1_0.xyz'], [[1,0,0], [0,1,0]])

    mergedGrid = pysurf.plot3d_interface.read_plot3d('merged.xyz',3)

    mergedGrid.remove_curves()

    pysurf.plot3d_interface.export_plot3d(mergedGrid, 'merged.xyz', saveNumpy=True)

    # Copy the numpy array containing coordinate info from the surface mesh
    # of the collar so we can access it later without recreating the mesh.
    subprocess.call(["cp merged.npy merged_"+str(i).zfill(2)+".npy"], shell=True)


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

    subprocess.call(["cp collar.cgns collar_master.cgns"], shell=True)

    # Set options for the pywarpustruct instance
    options = {
      'gridFile':'collar_master.cgns',
      'fileType':'CGNS',
      'specifiedSurfaces':None,
      'symmetrySurfaces':None,
      'symmetryPlanes':[],
      'aExp': 3.0,
      'bExp': 5.0,
      'LdefFact':100.0,
      'alpha':0.25,
      'errTol':0.0001,
      'evalMode':'fast',
      'useRotations':True,
      'zeroCornerRotations':True,
      'cornerAngle':30.0,
      'bucketSize':8,
    }

    # Create the mesh object using pywarpustruct
    mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

    # Get the initial coordinates in the order that the mesh object
    # uses them
    warp_coords = mesh.getSurfaceCoordinates()

    # Load the coordinates from the collar surface mesh we marched
    # earlier. Note that here in the first iteration the coordinate
    # points between this and the warp coords should be exactly the same.
    #
    # However, their order is not the same, so we must match the
    # coordinate points and save the indices so that we can reorder
    # our future surface marched points into the order that
    # pywarpustruct expects.
    march_coords = np.load('merged.npy')

    # To do this, we setup a KDTree using the pySurf generated
    # surface mesh coordinates.
    tree = cKDTree(march_coords)

    # We then query this tree with the warp coordinates to obtain
    # `index`, the mapping of the coordinate points.
    d, index = tree.query(warp_coords)

    return mesh, index

def run_pywarp_for_collar(mesh, index):

    print 'Using previously created volume mesh and warping it to the new surface.\n'
    coords = np.load('merged.npy')

    # Here we use the remapped coordinate points to pass in to
    # pywarpustruct.
    mesh.setSurfaceCoordinates(coords[index])

    # Actually warp the mesh and then write out the new volume mesh.
    mesh.warpMesh()
    mesh.writeGrid('collar.cgns')

    dXvWarpb = np.random.random(mesh.warp.griddata.warpmeshdof)
    mesh.warpDeriv(dXvWarpb, solverVec=False)
    # print 'here are the dXs!'
    # for i in mesh.getdXs():
    #     print i


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

    subprocess.call(["cgns_utils combine cube_vol.cgns rect_vol_temp.cgns collar.cgns half_full.cgns"], shell=True)

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
        'writeTecplotSurfaceSolution':True,
    }

    # Create solver
    CFDSolver = AEROSOLVER(options=aeroOptions, comm=comm, debug=True)

    # Here we just want to check flooding
    CFDSolver.setAeroProblem(ap)
    CFDSolver.writeSolution()

    subprocess.call(["cp fc_-001_surf.plt fc_surf_"+str(i).zfill(2)+".plt"], shell=True)

n = 2
for i in range(1):


    # These settings are known to be probably incorrect for valid overset
    # mesh generation, but we are using them primarily to check derivative
    # seed passing in and out of pySurf.
    extent = 0.0001
    trans = (i) * extent / (n - 1)
    rectTranslation = np.array([trans, trans, trans])

    extrude_cube_volume_mesh()
    extrude_rect_volume_mesh()

    # Here we actually run all the functions we just defined.
    march_surface_meshes()
    merge_surface_meshes()

    if i == 0:
        mesh, index = run_pyhyp_for_collar()
    else:
        run_pywarp_for_collar(mesh, index)

    create_OCart_mesh_and_merge()
    run_ADflow_to_check_connections()
