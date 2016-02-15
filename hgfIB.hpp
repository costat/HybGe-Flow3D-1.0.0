#include <vector>

void
immersedBoundary ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                   std::vector<int>& matJs, std::vector<double>& matVals );

void
immersedBoundarySingleComponent ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                                  std::vector<int>& matJs, std::vector<double>& matVals, \
                                  int component );

void BuildImmersedBoundary( FluidMesh& Mesh, double vf, int nObs );
