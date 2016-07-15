rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -c -cpp precision.f90

mpifort -c -cpp adtData.f90
mpifort -c -cpp gradD2Hexa.f90
mpifort -c -cpp hessD2Hexa.f90
mpifort -c -cpp minD2Hexa.f90
mpifort -c -cpp newtonStep.f90

mpifort -c -cpp adtUtils.f90

mpifort -c -cpp adtBuild.f90
mpifort -c -cpp adtLocalSearch.f90

mpifort -c -cpp adtSearch.f90

mpifort -c -cpp adtAPI.f90

mpifort *.o test.f90 -o test.exe

mpirun -np 1 ./test.exe
