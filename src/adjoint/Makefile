HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = adStack.o \
                 utilities_b.o \
                 utilities_d.o

FF90_OBJECTS_2 = intersection_b.o \
                 intersection_d.o \
                 adtProjections_b.o \
                 adtProjections_d.o \
                 curveUtils_d.o \
                 curveUtils_b.o

default: all

all:
	$(PMAKE) $(FF90_OBJECTS_1)
	$(PMAKE) $(FF90_OBJECTS_2)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in Intersection Derivatives... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
