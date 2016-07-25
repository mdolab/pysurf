include ${CONFIG_FILE}
include ${RULES_FILE}

FF90_OBJECTS_1 = adtData.o\
		 gradD2Hexa.o\
		 hessD2Hexa.o\
		 minD2Hexa.o\
		 newtonStep.o

FF90_OBJECTS_2 = adtUtils.o

FF90_OBJECTS_3 = adtBuild.o\
		 adtLocalSearch.o

FF90_OBJECTS_4 = adtSearch.o

FF90_OBJECTS_5 = adtAPI.o

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	$(PMAKE)  src3
	$(PMAKE)  src4
	$(PMAKE)  src5

clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)

src3:	$(FF90_OBJECTS_3)

src4:	$(FF90_OBJECTS_4)

src5:	$(FF90_OBJECTS_5)
