"""
Script to produce vertically translated wing-body meshes using pywarpustruct
after producing the initial collar mesh using pyHyp.
"""

import subprocess
import numpy as np
from mpi4py import MPI
from pywarpustruct import USMesh
from scipy.spatial import cKDTree
from time import time

# Set options for the collar pywarpustruct instance
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

# Set options for the wing pywarpustruct instance
wing_options = {
  'gridFile':'primary_meshes/wing_L1_hyp.cgns',
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

# Create a dictionary to store timing information
times = {}
times['intersections'] = []
times['mergeCollar'] = []
times['pyhyp'] = []
times['pywarp'] = []
times['mergeMeshes'] = []
times['ADflow'] = []

# Function to print time formatted info
def print_line_time(name, time_list, string=''):
    sumlist = np.sum(time_list)
    print('   %-14s : %10.7f : %10.7f' % (name, sumlist, sumlist / len(time_list)), string)

if __name__ == "__main__":

    nStates = 2

    wingTranslation = np.zeros((nStates,3))
    wingTranslation[:, 2] = np.linspace(0., 3., nStates)

    for i, wingT in enumerate(wingTranslation):
        # print "Extracting curves.\n"
        # subprocess.call(["python script01_curveExtraction.py"], shell=True)
        #
        # print "Fixing trailing edge.\n"
        # subprocess.call(["python script02_fix_trailing_edge.py"], shell=True)

        print("Replacing values in scripts with desired translation.\n")
        python_string_replacement = "sed -i -e 's/wingTranslation =.*/wingTranslation = [{}, {}, {}]/g' script03_example_crm_fullInt.py".format(wingT[0], wingT[1], wingT[2])
        subprocess.call([python_string_replacement], shell=True)

        cgns_string_replacement = "sed -i -e 's/cgns_utils coarsen crm_wb.cgns .*/cgns_utils coarsen crm_wb.cgns crm_wb_L2_" + str(i).zfill(2) + ".cgns/g' script06_alternative.sh"
        subprocess.call([cgns_string_replacement], shell=True)

        t1 = time()

        print("Computing intersections and marching collar mesh.\n")
        subprocess.call(["python script03_example_crm_fullInt.py"], shell=True)

        t2 = time()
        times['intersections'].append(t2 - t1)

        print("Merging collar meshes.\n")
        subprocess.call(["python script04_mergeCollar.py"], shell=True)

        t3 = time()
        times['mergeCollar'].append(t3 - t2)

        # Copy the numpy array containing coordinate info from the surface mesh
        # of the collar so we can access it later without recreating the mesh.
        subprocess.call(["cp merged.npy merged_"+str(i).zfill(2)+".npy"], shell=True)

        # In the first iteration, save the first mesh that is produced by
        # pyHyp so we can later warp it
        if i == 0:

            print("Extruding collar mesh using pyHyp.\n")
            subprocess.call(["python script05_runPyhyp.py"], shell=True)

            t4 = time()
            times['pyhyp'].append(t4 - t3)

            subprocess.call(["cp collar.cgns collar_master.cgns"], shell=True)

            # Create the mesh object using pywarpustruct
            collar_mesh = USMesh(options=options, comm=MPI.COMM_WORLD)

            subprocess.call(["cp primary_meshes/wing_L1_hyp.cgns wing_L1_temp.cgns"], shell=True)

            # Create the wing mesh object using pywarpustruct
            wing_mesh = USMesh(options=wing_options, comm=MPI.COMM_WORLD)

            # Get the wing coordinates from the primary mesh
            wing_coords = wing_mesh.getSurfaceCoordinates()

            # Get the initial coordinates in the order that the mesh object
            # uses them
            warp_coords = collar_mesh.getSurfaceCoordinates()

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

            t5 = time()
            times['pywarp'].append(t5 - t4)

        else:

            print('Using previously created volume mesh and warping it to the new surface.\n')
            coords = np.load('merged.npy')

            # Here we use the remapped coordinate points to pass in to
            # pywarpustruct.
            collar_mesh.setSurfaceCoordinates(coords[index])

            # Actually warp the mesh and then write out the new volume mesh.
            collar_mesh.warpMesh()
            collar_mesh.writeGrid('collar.cgns')

            # Translate the initial wing coords based on the wingT values
            new_wing_coords = wing_coords + wingT / .0254

            # Set the new wing coords, warp the mesh, and save the grid
            wing_mesh.setSurfaceCoordinates(new_wing_coords)
            wing_mesh.warpMesh()
            wing_mesh.writeGrid('wing_L1_temp.cgns')

            t5 = time()
            times['pywarp'].append(t5 - t3)

        print("Merging meshes into single .cgns file.\n")
        subprocess.call(["sh script06_alternative.sh"], shell=True)

        t6 = time()
        times['mergeMeshes'].append(t6 - t5)

        print("Renaming and copying meshes.\n")
        subprocess.call(["cp crm_wb_L1.cgns grids/crm_wb_L1_"+str(i).zfill(2)+".cgns"], shell=True)
        subprocess.call(["cp crm_wb_L2.cgns grids/crm_wb_L2_"+str(i).zfill(2)+".cgns"], shell=True)

        print("Running ADflow using the combined meshes.\n")
        subprocess.call(["cd ADflow; sh run_check.sh"], shell=True)
        subprocess.call(["cd ADflow; cp fc_-001_surf.plt fc_surf_"+str(i).zfill(2)+".plt"], shell=True)

        t7 = time()
        times['ADflow'].append(t7 - t6)

    print('           Timings summary (sec)      ')
    print('   ========================================  ')
    print('                       Total       Average')
    print('                                          ')

    for key, value in times.iteritems():
        print_line_time(key, value)
    print()
