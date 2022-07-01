HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = cgnsGrid.o

FF90_OBJECTS_2 = cgnsInterface.o

FF90_OBJECTS_3 = cgnsAPI.o

FILES_TO_COMPLEXIFY = cgnsGrid.F90 \
                      cgnsInterface.F90 \
                      cgnsAPI.F90

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	$(PMAKE) $(FF90_OBJECTS_3)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in CGNSInterface... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/CGNSInterface/$$file -o ! -e ../../src_cs/CGNSInterface/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/CGNSInterface/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/CGNSInterface
