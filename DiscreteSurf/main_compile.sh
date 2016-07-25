cd src/common
bash compile_and_run.sh

cd ../ADT
bash compile_and_run.sh

cd ../CGNSInterface
bash compile_and_run.sh

cd ../../lib
ar -rv libdiscsurf.a ../obj/*.o

cd ../src/examples
bash compile_and_run.sh
./test.exe
./testADT.exe
./test2.exe
cd ../..
