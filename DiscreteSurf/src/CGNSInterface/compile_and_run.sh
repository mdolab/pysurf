rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp -g -fbounds-check cgnsGrid.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src CGNSinterface.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

#mpifort -I../../mod -c -cpp -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src readCGNS.F90

#mpifort *.o -L$HOME/packages/cgnslib_3.2.1/src -lcgns -g -fbounds-check test.f90 -o test.exe

#mpirun -np 1 ./test.exe

#mv *.o ../../obj/.
#mv *.mod ../../mod/.
