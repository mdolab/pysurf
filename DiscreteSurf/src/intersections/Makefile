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

FF90_OBJECTS_1 = Intersection.o \

FF90_OBJECTS_2 = intersectionAPI.o \
                 test.o

default: all

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in Intersections... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)
