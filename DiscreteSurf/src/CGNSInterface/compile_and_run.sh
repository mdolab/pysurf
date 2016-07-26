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

mpifort -I../../mod -c -fPIC CGNSapi.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

f2py CGNSapi.F90 -m CGNSapi -h CGNSapi.pyf --overwrite-signature

f2py -c -L../../lib -ldiscsurf -L$HOME/packages/cgnslib_3.2.1/src -lcgns -L/home/$USER/packages/petsc-3.6.1/real-debug/lib -lmpi_usempi -lmpi_mpifh -lmpi -I../../mod ../../obj/CGNSapi.o CGNSapi.pyf

python importTest.py

mv *.so ../../python/.
