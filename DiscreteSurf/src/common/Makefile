include ${CONFIG_FILE}
include ${RULES_FILE}

FF90_OBJECTS_1 = precision.o

FF90_OBJECTS_2 = constants.o

FF90_OBJECTS_3 = cgnsGrid.o

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	$(PMAKE)  src3

clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)

src3:	$(FF90_OBJECTS_3)
