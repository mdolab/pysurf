#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: Ney Secco and John Jasa                               *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 02-14-2017                                      *
#      * Last modified: 02-14-2017                                      *
#      *                                                                *
#      ******************************************************************

HOME_DIR = ..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = precision.o

FF90_OBJECTS_2 = structCGNS_API.o


default: all

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in CGNSInterface... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)
