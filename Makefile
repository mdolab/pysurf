#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: Ney Secco and John Jasa                               *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 07-27-2016                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = geometryEngines/TSurf \
		meshTools/hypsurf \
	        utilities/CGNSinterface \

default:
# Check if the config.mk file is in the config dir.
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done

clean:
	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
