#!/bin/bash

# Need one argument at least
if [ -z $1 ]; then
    echo "$0 needs to be called with the compiler name as an argument"; exit 1
fi

# Initialize variables
FF90=$1
FFLAGS=""

# Allow argument mismatch for gfortran >= 10

# NOTE: The reason for adding this flag is because Tapenade 3.16's ADFirstAidKIt is not compliant with GCC 10+.
# This should be removed once Tapenade fixes this or removes the ADFirstAidKIt.

fc=$("$FF90" --version 2>&1 | grep -i 'gnu')
if [ ! -z "$fc" ]; then
    # Get the GNU compiler version
    version=$("$FF90" -v 2>&1 | grep 'gcc version' | cut -d' ' -f3)
    if [ ! -z "$version" ]; then
      # Get the major version
      version_major=`echo $version | cut -f1 -d.`
    fi

    if [ $version_major -ge 10 ]; then
        FFLAGS="-fallow-argument-mismatch"
    fi
fi

# Print at end to add to the makefile
echo "$FFLAGS"
