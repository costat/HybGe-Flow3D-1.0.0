#include <vector>

#include <paralution.hpp>

#include "hgfMesh.hpp"

void
StokesSolveDirect( const FluidMesh& Mesh, double visc, int direction, \
                   std::vector<double>& Solution, double tolAbs, double tolRel, \
                   int maxIt, int nThreads, int prec );
void
StokesSolveRich( const FluidMesh& Mesh, double visc, int direction, \
                 std::vector<double>& Solution, double tolAbs, double tolRel, \
                 int maxIt, int nThreads, int prec );
void
initPressure( const FluidMesh& Mesh, std::vector<double>& Solution );
void
setForceRich( const FluidMesh& Mesh, const std::vector<double>& Solution, std::vector<double>& b, \
              const std::vector<double>& force, int component );
void
updatePressureRich( const FluidMesh& Mesh, std::vector<double>& Solution, \
                    const paralution::LocalVector<double>& solU, const paralution::LocalVector<double>& solV, \
                    double& residual );
