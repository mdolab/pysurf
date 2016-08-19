#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: Ney Secco and John Jasa                               *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 07-27-2016                                      *
#      *                                                                *
#      ******************************************************************
#
HOME_DIR = ..
COMMON_FILE = $(HOME_DIR)/common.mk
include ${COMMON_FILE}

default: all

all: 	libhypsurf

clean:
	@echo "        Making clean in lib ... "
	rm -f lib/*

libhypsurf:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) libhypsurf.a $(OBJDIR)/*.o
	@echo
	@echo "        library libhypsurf.a created."
	@echo
