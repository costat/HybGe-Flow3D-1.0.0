#include <vector>

#include "hgfMesh.hpp"

void
StokesSolve( const FluidMesh& Mesh, double visc, int direction, \
             std::vector<double>& Solution, double tolAbs, double tolRel, \
             int maxIt, int nThreads, int prec );

