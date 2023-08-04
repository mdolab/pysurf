# ----------------------------------------------------------------------
# Config file for ifort
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpiifort
CC   = mpiicc

# ----------- CGNS ------------------
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# ------- Define complexify inlcude and linker flags -------------------------
COMPLEXIFY_INCLUDE_FLAGS=-I$(COMPLEXIFY_DIR)/include
COMPLEXIFY_LINKER_FLAGS=-L$(COMPLEXIFY_DIR)/lib -lcomplexify

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -r8 -O2 -g
FF90_FLAGS = ${FF77_FLAGS} -std08
C_FLAGS    = -fPIC -O2

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER = $(FF90)
LINKER_FLAGS = -nofor-main

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py
