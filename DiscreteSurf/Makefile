# Makefile for DiscSurf
# Ney Secco 07-22-2016
# 

# Specify calling options
.PHONY: clean help

# Specify config and rules file
CONFIG_FILE = config.mk
RULES_FILE = common.mk

# Folder with source file where we should run make
SRC_DIR = src/modules \
          src/ADT \
          src/CGNSinterface

# Important directories
MOD_DIR = $(CURDIR)/mod
OBJ_DIR = $(CURDIR)/obj

# Add config and rules file
include ${CONFIG_FILE}
include ${RULES_FILE}

default:
# Run make in each subdirectory
	@for subdir in $(SRC_DIR) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			export MOD_DIR; \
			export OBJ_DIR; \
			export CONFIG_FILE; \
			export RULES_FILE; \
			$(MAKE) all -C $$subdir; \
			(cd $$subdir && mv *.o $(OBJ_DIR)/.); \
			(cd $$subdir && mv *.mod $(MOD_DIR)/.); \
		done

