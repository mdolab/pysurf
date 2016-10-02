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

FF90_OBJECTS_1 = adBuffer.o \
		 adStack.o \
		 outputReverse/hypsurfmain_b.o \
		 outputForward/hypsurfmain_d.o

default: all

all:
# TODO: remove these calls to tapenade here
# Make the forward mode derivative files
	make -f Makefile_tapenade

# Make the reverse mode derivative files
	make -f Makefile_tapenade ad_reverse

	$(PMAKE)  src1
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in HypsurfAPI Derivatives... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)
