nvcc -v ../src/driver.cpp ../src/hgfM.cpp ../src/hgfArrays.cpp ../src/hgfPPM.cpp ../src/hgfMeshCu.cu ../src/hgfBC.cpp ../src/hgfIB.cpp ../src/hgfPoreNetworkM.cpp ../src/hgfStokesM.cpp -O3 -o hgfTest -Xcompiler -fopenmp -lcuda -DCUDA_BUILD=1 -I../src/
