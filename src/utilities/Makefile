#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: Ney Secco and John Jasa                               *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 02-01-2017                                      *
#      *                                                                *
#      ******************************************************************

HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_0 = complexify.o \

FF90_OBJECTS_1 = Utilities.o \

FF90_OBJECTS_2 = utilitiesAPI.o \

FILES_TO_COMPLEXIFY =	Utilities.F90\
			utilitiesAPI.F90\

FILES_TO_COPY = complexify.f90\

default: all

all:
	$(PMAKE)  src0
	$(PMAKE)  src1
	$(PMAKE)  src2
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in Utilities... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/utilities/$$file -o ! -e ../../src_cs/utilities/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/utilities/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/utilities/$$file -o ! -e ../../src_cs/utilities/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/utilities/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/utilities


src0: $(FF90_OBJECTS_0)
src1:	$(FF90_OBJECTS_1)
src2:	$(FF90_OBJECTS_2)