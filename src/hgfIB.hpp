// hgfIB.hpp
#ifndef HGFIB_H
#define HGFIB_H

#include <vector>

void
immersedBoundary ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                   std::vector<int>& matJs, std::vector<double>& matVals );

void
immersedBoundarySingleComponent ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                                  std::vector<int>& matJs, std::vector<double>& matVals, \
                                  int component );

int BuildImmersedBoundary( FluidMesh& Mesh, double vf, int nObs );

#endif
