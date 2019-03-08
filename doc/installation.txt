### Build instructions

Before you build anything, you need to adjust the config files for the following packages:

- geometryEngines/TSurf
- meshTools/hypsurf
- utilities/CGNSinterface

Navigate to each of these folder and make a copy of the appropriate config file.
For instance:
$ cd geometryEngines/TSurf
$ cp config/defaults/config_LINUX_GFORTRAN.mk config/config.mk

Once you adjust all config files, go back to the root repo in pySurf and run:
$ make

This should build all relevant packages.
