rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

#mpifort -I../../mod -c -cpp -g -fbounds-check test.F90
mpifort -I../../mod -c -cpp -g -fbounds-check testADT.F90
#mpifort -I../../mod -c -cpp -g -fbounds-check test2.F90

#mpifort test.o -L../../lib -ldiscretesurf -L$HOME/packages/cgnslib_3.2.1/src -lcgns -o test.exe

mpifort testADT.o -L../../lib -ldiscretesurf -L$HOME/packages/cgnslib_3.2.1/src -lcgns -o testADT.exe

#mpifort test2.o -L../../lib -ldiscretesurf -L$HOME/packages/cgnslib_3.2.1/src -lcgns -o test2.exe

#mv *.o ../../obj/.
#mv *.mod ../../mod/.
