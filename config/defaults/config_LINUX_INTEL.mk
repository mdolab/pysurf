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

# ------- Define Compiler Flags ----------------------------------------
FF90_GEN_FLAGS = -DHAS_ISNAN -fPIC -r8 -O2 -g
CC_GEN_FLAGS   = -DHAS_ISNAN -O -fPIC

FF90_OPT_FLAGS   = -DHAS_ISNAN -fPIC -r8 -O2 -g
CC_OPT_FLAGS     = -DHAS_ISNAN -O -fPIC

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER = $(FF90)
LINKER_FLAGS = -nofor_main

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py
