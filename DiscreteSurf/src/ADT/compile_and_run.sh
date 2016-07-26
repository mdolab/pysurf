rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp -fPIC adtData.F90
mpifort -I../../mod -c -cpp -fPIC gradD2Hexa.F90
mpifort -I../../mod -c -cpp -fPIC hessD2Hexa.F90
mpifort -I../../mod -c -cpp -fPIC minD2Hexa.F90
mpifort -I../../mod -c -cpp -fPIC newtonStep.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC adtUtils.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC adtBuild.F90
mpifort -I../../mod -c -cpp -fPIC adtLocalSearch.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC adtSearch.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp -fPIC adtAPI.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

f2py -c -L../../lib -ldiscsurf -L/home/$USER/packages/petsc-3.6.1/real-debug/lib -lmpi_usempi -lmpi_mpifh -lmpi -I../../mod ../../obj/adtAPI.o adtAPI.pyf

python importTest.py

mv *.so ../../python/.
