HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = precision.o

FF90_OBJECTS_2 = communication.o \
                 constants.o

FILES_TO_COMPLEXIFY = precision.F90 \
                      communication.F90 \
                      constants.F90

FILES_TO_COPY =

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in common... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/common/$$file -o ! -e ../../src_cs/common/$$file ]; then \
			echo complexifying file $$file; \
			complexify $$file; \
			mv c_$$file ../../src_cs/common/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/common/$$file -o ! -e ../../src_cs/common/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/common/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/common
