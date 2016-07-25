rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp -fPIC -g -fbounds-check cgnsGrid.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src CGNSinterface.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

# mpifort -I../../mod -c -fPIC CGNSinter.f90

f2py CGNSinterface.F90 -m CGNSinterface -h CGNSinterface.pyf --overwrite-signature

f2py -lgfortran -c -I../../mod CGNSinterface.o CGNSinterface.pyf
