# General variables and rules to make objects

# Assemble general flags
FF90_ALL_FLAGS = -I$(MOD_DIR) $(CGNS_INCLUDE_FLAGS) $(FF90_FLAGS)

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo
