rm *.o
rm *.exe
reset

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp -fPIC -g -fbounds-check cgnsGrid.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC -g -fbounds-check -I$HOME/packages/cgnslib_3.2.1/src CGNSinterface.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -fPIC CGNSf2py.f90

mv *.o ../../obj/.
mv *.mod ../../mod/.

f2py CGNSf2py.f90 -m CGNSf2py -h CGNSf2py.pyf --overwrite-signature

f2py -c -L../../lib -ldiscsurf -L$HOME/packages/cgnslib_3.2.1/src -lcgns -L/home/$USER/packages/petsc-3.6.1/real-debug/lib -lmpi_usempi -lmpi_mpifh -lmpi -I../../mod ../../obj/CGNSf2py.o CGNSf2py.pyf

python importTest.py
