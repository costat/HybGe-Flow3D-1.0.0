#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <amgx_c.h>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

// solves the stokes system directly as a single linear system
void
StokesSolveDirect( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par )
{

  // delcarations
  std::vector<int> matIs, matJs, rowPTR;;
  std::vector<double> matVals, force;
  matIs.reserve(Mesh.maxNNZ);
  matJs.reserve(Mesh.maxNNZ);
  matVals.reserve(Mesh.maxNNZ);
  force.resize(Mesh.dofTotal);

  // interior cells
  StokesArray( Mesh, Par.visc, matIs, matJs, matVals );

  // boundary conditions
  AxisFlowDrive( Mesh, matIs, matJs, matVals, force, Par.visc, Par.direction );

  // immersed boundary
  immersedBoundary( Mesh, matIs, matJs, matVals );

  // build the rowPTR vector for CSR rep of the array. first we sort the COO vecs
  sortCOO( matIs, matJs, matVals );
  buildCSR( matIs, matJs, matVals, rowPTR );

  //---- CREATE AMGX OBJECTS ----//
  // AMGX config object
  AMGX_config_handle config;
  AMGX_config_create_from_file( &config, "configs/PBICGSTAB.json");

  // AMGX resources object
  AMGX_resources_handle resources;
  AMGX_resources_create_simple(&resources, config);

  // AMGX solver object
  AMGX_solver_handle solver;
  AMGX_solver_create(&solver, resources, AMGX_mode_dDDI, config);

  // AMGX matrix object
  AMGX_matrix_handle matrix;
  AMGX_matrix_create(&matrix, resources, AMGX_mode_dDDI);

  // AMGX vector objects
  AMGX_vector_handle x, b;
  AMGX_vector_create(&x, resources, AMGX_mode_dDDI);
  AMGX_vector_create(&b, resources, AMGX_mode_dDDI);

  // Build matrix, b, x
  AMGX_matrix_upload_all(matrix, Mesh.dofTotal, matVals.size(), 1, 1, rowPTR.data(), matJs.data(), matVals.data(), 0);
  AMGX_vector_upload(b, Mesh.dofTotal, 1, force.data());
  AMGX_vector_set_zero(x, Mesh.dofTotal, 1);

  // Setup solver from matrix
  AMGX_solver_setup(solver, matrix);

  // solve the system with 0 initial guess
  AMGX_solver_solve(solver, b, x);

  // copy solution to std::vector reference
  AMGX_vector_download(x, Solution.data());

  // cleanup
  AMGX_solver_destroy(solver);
  AMGX_matrix_destroy(matrix);
  AMGX_vector_destroy(x);
  AMGX_vector_destroy(b);
  AMGX_resources_destroy(resources);
  AMGX_config_destroy(config);

}

void SolverInit( void )
{
  AMGX_initialize();
  AMGX_initialize_plugins();
}

void SolverFinalize( void )
{
  AMGX_finalize_plugins();
  AMGX_finalize();
}
