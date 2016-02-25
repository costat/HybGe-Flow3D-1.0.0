// hgfStokesM.hpp
#ifndef HGFSTOKES_H
#define HGFSTOKES_H

#include <vector>

// hgf includes
#include "hgfMeshCu.hpp"

void
StokesSolveDirect( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par );
void
SolverInit( void );
void
SolverFinalize( void );

#endif
