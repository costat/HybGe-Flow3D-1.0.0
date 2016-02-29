// hgfStokesM.hpp
#ifndef HGFSTOKES_H
#define HGFSTOKES_H

#include <vector>
#include <paralution.hpp>

// hgf includes
#include "hgfMeshCu.hpp"

void
StokesSolveDirect( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par );
void
InitPressure( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par );
void
SetForceUZCG( const FluidMesh& Mesh, const std::vector<double>& Solution, LocalVector<double>& b, \
                                     const std::vector<double>& force, int component );
void
StokesSolveUZCG( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par );
void
SolverInit( void );
void
SolverFinalize( void );

#endif
