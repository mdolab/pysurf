rm *.o
rm *.exe

# The -cpp flag indicates that we should use the preprocessor to
# parse the #ifdef statements

mpifort -I../../mod -c -cpp adtData.F90
mpifort -I../../mod -c -cpp gradD2Hexa.F90
mpifort -I../../mod -c -cpp hessD2Hexa.F90
mpifort -I../../mod -c -cpp minD2Hexa.F90
mpifort -I../../mod -c -cpp newtonStep.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp adtUtils.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp adtBuild.F90
mpifort -I../../mod -c -cpp adtLocalSearch.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp adtSearch.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.

mpifort -I../../mod -c -cpp adtAPI.F90

mv *.o ../../obj/.
mv *.mod ../../mod/.
