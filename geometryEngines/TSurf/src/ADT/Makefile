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

FF90_OBJECTS_1 = adtData.o\
	         gradD2Hexa.o\
	         hessD2Hexa.o\
	         minD2Hexa.o\
	         newtonStep.o\

FF90_OBJECTS_2 = adtUtils.o \
                 adtProjections.o

FF90_OBJECTS_3 = adtBuild.o\
	         adtLocalSearch.o

FF90_OBJECTS_4 = adtSearch.o

FF90_OBJECTS_5 = adtAPI.o


default: all

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	$(PMAKE)  src3
	$(PMAKE)  src4
	$(PMAKE)  src5
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in ADT... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)

src3:	$(FF90_OBJECTS_3)

src4:	$(FF90_OBJECTS_4)

src5:	$(FF90_OBJECTS_5)
