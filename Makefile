SUBDIR_SRC    = src/common \
                src/adjoint \
                src/utilities \
                src/CGNSInterface \
                src/ADT \
                src/curveSearch \
                src/intersections \

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config_LINUX_GFORTRAN.mk config/config.mk"; \
	echo " ";\
	echo "Modify this config file as required. With the config file "; \
	echo "specified, rerun 'make' and the build will start."; \
	else \
		make discretesurf || exit 1; \
		make -f Makefile_CS discretesurf || exit 1; \
	fi;

clean:
	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*
	(cd pysurf && rm *.so) || exit 1;

	make -f Makefile_CS clean

discretesurf:
	mkdir -p obj
	mkdir -p mod
	ln -sf config/config.mk config.mk
	ln -sf common_real.mk common.mk
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)
	(cd src/f2py && make)
