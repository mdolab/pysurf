"""
This script will extract sharp corners (such as trailing edge) from the crm geometry
and save them in tecplot format.

Ney Secco 2016-11
"""

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import hashlib
import unittest
import os
import pickle


# TESTING FUNCTION

os.system("rm curve*")

wing = pysurf.TSurfGeometry("../inputs/rect_uns.cgns", MPI.COMM_WORLD)

wing.extract_curves()

for (curveID, curve) in enumerate(wing.curves.itervalues()):
    curve.export_tecplot("extracted_curve_%03d" % curveID)
