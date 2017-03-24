If you want to regenerate the f2py signature file use:

$ bash getf2pySignature.sh

But remeber to go to ../python/f2py/cgnsAPI.pyf and remove all allocatable variable declarations since this might have issues with the Intel compiler.
