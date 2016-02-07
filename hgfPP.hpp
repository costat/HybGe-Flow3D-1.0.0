#include <vector>
#include <string>

void
computeKTensorL ( const FluidMesh& Mesh, \
                  const std::vector<double>& xSolution, \
                  const std::vector<double>& ySolution, \
                  const std::vector<double>& zSolution );

void
computeAveragesX ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution,
                   double& V, double& G, int print );

void
computeAveragesY ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print );

void
computeAveragesZ ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print );

void
computeKConstantDrive ( const FluidMesh& Mesh, \
                        const std::vector<double>& Solution, \
                        int direction );

void
writeSolutionTP ( const FluidMesh& Mesh, const std::vector<double>& sol, \
                  std::string& outName );
