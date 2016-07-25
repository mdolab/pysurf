include ${CONFIG_FILE}
include ${RULES_FILE}

FF90_OBJECTS =	readUnstructuredCGNS.o\
		readCGNS.o

all:
	$(PMAKE)  src

clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src:	$(FF90_OBJECTS)
