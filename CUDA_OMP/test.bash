nvcc -std=c++11 ../src/driver.cpp ../src/hgfM.cpp ../src/hgfArrays.cpp ../src/hgfPPM.cpp ../src/hgfMeshCu.cu ../src/hgfBC.cpp ../src/hgfIB.cpp ../src/hgfPoreNetworkM.cpp ../src/hgfStokesM.cpp -O3 -o hgfTest -lcuda -lcusparse -lcublas -Xcompiler -fopenmp -lmagma -lmagma_sparse -lblas -llapack -DCUDA_BUILD=1 -I../src/