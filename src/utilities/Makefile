HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = utilities.o

FF90_OBJECTS_2 = utilitiesAPI.o

FILES_TO_COMPLEXIFY = utilities.F90 \
                      utilitiesAPI.F90

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in utilities... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/utilities/$$file -o ! -e ../../src_cs/utilities/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/utilities/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/utilities
