#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: Ney Secco and John Jasa                               *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 07-27-2016                                      *
#      *                                                                *
#      ******************************************************************

HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = cgnsGrid.o

FF90_OBJECTS_2 = cgnsInterface.o

FF90_OBJECTS_3 = cgnsAPI.o


default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	$(PMAKE) $(FF90_OBJECTS_3)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in CGNSInterface... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
