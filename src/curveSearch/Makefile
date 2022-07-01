HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = curveUtils.o

FF90_OBJECTS_2 = curveSearchAPI.o

FILES_TO_COMPLEXIFY = curveUtils.F90 \
                      curveSearchAPI.F90

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in curveSearch... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/curveSearch/$$file -o ! -e ../../src_cs/curveSearch/$$file ]; then \
			echo complexifying file $$file; \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/curveSearch/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/curveSearch
