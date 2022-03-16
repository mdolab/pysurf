#      ******************************************************************
#      *                                                                *
#      * Description: Rules to make the objects. These are the general  *
#      *              rules. If in a subdirectory different rules must  *
#      *              be used, this file should not be included.        *
#      *                                                                *
#      ******************************************************************

%.o: %.F90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

%.o: %.f90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o: %.f
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

%.o: %.F
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

%.o: %.c
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
