// hgfStokesM.hpp
#ifndef HGFSTOKES_H
#define HGFSTOKES_H

#include <vector>

// hgf includes
#include "hgfMeshCu.cuh"

void
StokesSolveDirect( const FluidMesh& Mesh, double visc, int direction, \
                   std::vector<double>& Solution, double tolAbs, double tolRel, \
                   int maxIt, int nThreads, int prec );
void
StokesSolveRich( const FluidMesh& Mesh, double visc, int direction, \
                 std::vector<double>& Solution, double tolAbs, double tolRel, \
                 int maxIt, int nThreads, int prec, double relax );
void
initPressure( const FluidMesh& Mesh, std::vector<double>& Solution, int direction );
void
setForceRich( const FluidMesh& Mesh, const std::vector<double>& Solution, std::vector<double>& b, \
              const std::vector<double>& force, int component );
void
updatePressureRich( const FluidMesh& Mesh, std::vector<double>& Solution, \
                    const std::vector<double>& solU, const std::vector<double>& solV, \
                    const std::vector<double>& solW, double& residual, double relax );

#endif
