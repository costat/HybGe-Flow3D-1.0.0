#include <vector>
#include <paralution.hpp>

void
computeKTensorL ( const FluidMesh& Mesh, \
                  const paralution::LocalVector<double>& xSolution, \
                  const paralution::LocalVector<double>& ySolution, \
                  const paralution::LocalVector<double>& zSolution );

void
computeAveragesX ( const FluidMesh& Mesh, \
                   const paralution::LocalVector<double>& Solution,
                   double& V, double& G, int print );

void
computeAveragesY ( const FluidMesh& Mesh, \
                   const paralution::LocalVector<double>& Solution, \
                   double& V, double& G, int print );

void
computeAveragesZ ( const FluidMesh& Mesh, \
                   const paralution::LocalVector<double>& Solution, \
                   double& V, double& G, int print );

void
computeKConstantDrive ( const FluidMesh& Mesh, \
                        const paralution::LocalVector<double>& Solution, \
                        int direction );

void
writeSolutionL ( const FluidMesh& Mesh, const paralution::LocalVector<double>& sol, \
                 std::string& outName );

