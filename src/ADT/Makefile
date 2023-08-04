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
                 quadProjDiff.o

FF90_OBJECTS_2 = adtUtils.o \
                 adtProjections.o

FF90_OBJECTS_3 = adtBuild.o\
                 adtLocalSearch.o

FF90_OBJECTS_4 = adtSearch.o

FF90_OBJECTS_5 = adtAPI.o

FILES_TO_COMPLEXIFY = adtData.F90 \
                      gradD2Hexa.F90 \
                      hessD2Hexa.F90 \
                      minD2Hexa.F90 \
                      newtonStep.F90 \
                      quadProjDiff.F90 \
                      adtUtils.F90 \
                      adtProjections.F90 \
                      adtBuild.F90 \
                      adtLocalSearch.F90 \
                      adtSearch.F90 \
                      adtAPI.F90

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	$(PMAKE) $(FF90_OBJECTS_3)
	$(PMAKE) $(FF90_OBJECTS_4)
	$(PMAKE) $(FF90_OBJECTS_5)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in ADT... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/ADT/$$file -o ! -e ../../src_cs/ADT/$$file ]; then \
			echo complexifying file $$file; \
			complexify $$file; \
			mv c_$$file ../../src_cs/ADT/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/ADT
