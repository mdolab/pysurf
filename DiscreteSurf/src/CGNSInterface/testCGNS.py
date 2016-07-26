# IMPORTS
import CGNSf2py
from mpi4py import MPI

# Set CGNS file
fileName = 'cube2.cgns'

# Set communicator
comm = MPI.COMM_WORLD

print CGNSf2py.cgnsf2py.__doc__

# Read CGNS file
CGNSf2py.cgnsf2py.cgnsf2py_routine(fileName,0)

print CGNSf2py.cgnsf2py.__doc__

# Print output
print CGNSf2py.cgnsf2py.quadsconn
