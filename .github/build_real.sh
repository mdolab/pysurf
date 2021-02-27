#!/bin/bash
set -e

cd geometryEngines/TSurf
cp $CONFIG_FILE config/config.mk
cd ../../meshTools/hypsurf
cp $CONFIG_FILE config/config.mk
cd ../../utilities/CGNSinterface
cp $CONFIG_FILE config/config.mk

cd ../..
make
