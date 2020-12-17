import pysurf

'''
This script takes the two blocks that define the collar mesh
and join them so that we can run pyHyp
'''

pysurf.plot3d_interface.merge_plot3d(['wing.xyz', 'body.xyz'], [[1,0,0], [0,1,0]])

mergedGrid = pysurf.plot3d_interface.read_plot3d('merged.xyz',3)

mergedGrid.remove_curves()

pysurf.plot3d_interface.export_plot3d(mergedGrid, 'merged.xyz', saveNumpy=True)
