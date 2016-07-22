rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -c -cpp -g -fbounds-check precision.f90
mpifort -c -cpp -g -fbounds-check communication.f90
mpifort -c -cpp -g -fbounds-check constants.F90
mpifort -c -cpp -g -fbounds-check cgnsGrid.F90

mpifort -c -cpp -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src readUnstructuredCGNS.F90
mpifort -c -cpp -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src readCGNS.F90

mpifort *.o -L$HOME/packages/cgnslib_3.2.1/src -lcgns -g -fbounds-check test.f90 -o test.exe

mpirun -np 1 ./test.exe
