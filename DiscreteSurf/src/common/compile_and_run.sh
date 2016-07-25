rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp -g -fbounds-check precision.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -g -fbounds-check communication.F90
mpifort -I../../mod -c -cpp -g -fbounds-check constants.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.
