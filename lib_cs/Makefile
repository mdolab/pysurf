HOME_DIR = ..
COMMON_FILE = $(HOME_DIR)/common.mk
include ${COMMON_FILE}

default: all

all: 	libdiscretesurf

clean:
	@echo "        Making clean in lib ... "
	rm -f lib/*

libdiscretesurf:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) libdiscretesurf_cs.a $(OBJDIR)/*.o
	@echo
	@echo "        library libdiscretesurf_cs.a created."
	@echo
