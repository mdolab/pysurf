#!/bin/bash
# This file will download the input files for the pySurf tests and extract them to the right place.

DIR=$(dirname $0)
TAR="input_files.tar.gz"
wget -O $DIR/$TAR https://websites.umich.edu/~mdolaboratory/repo_files/pySurf/pysurf_input_files.tar.gz
tar -xzf $DIR/$TAR -C $DIR/../
rm $DIR/$TAR
