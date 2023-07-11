MODDIR = $(HOME_DIR)/mod_cs
OBJDIR = $(HOME_DIR)/obj_cs
LIBDIR = $(HOME_DIR)/lib_cs

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

COMPILERS = $(HOME_DIR)/config.mk
ifneq ($(MAKECMDGOALS),clean)
include ${COMPILERS}
endif

#      ******************************************************************
#      *                                                                *
#      * Redefine .SUFFIXES to be sure all the desired ones are         *
#      * included.                                                      *
#      *                                                                *
#      ******************************************************************

.SUFFIXES: .o .f .F .f90 .F90

#      ******************************************************************
#      *                                                                *
#      * Arguments of make clean.                                       *
#      *                                                                *
#      ******************************************************************

MAKE_CLEAN_ARGUMENTS = *~ *.o *.mod *.il *.stb c_* *.a

#      ******************************************************************
#      *                                                                *
#      * Compiler flags to compile the sources.                         *
#      * The macro's ADDITIONAL_FF90_FLAGS and ADDITIONAL_CC_ALL_FLAGS  *
#      * make it possible that every subdirectory adds its own specific *
#      * compiler flags, if necessary.                                  *
#      *                                                                *
#      ******************************************************************

FF77_ALL_FLAGS = -I$(MODDIR) \
		$(CGNS_INCLUDE_FLAGS) \
		$(FF77_FLAGS)

FF90_ALL_FLAGS = -I$(MODDIR) \
		$(CGNS_INCLUDE_FLAGS) \
		$(COMPLEXIFY_INCLUDE_FLAGS) \
		$(FF90_FLAGS) \
		-DUSE_COMPLEX

CC_ALL_FLAGS = -I$(MODDIR) \
		$(CGNS_INCLUDE_FLAGS) \
		$(C_FLAGS)
