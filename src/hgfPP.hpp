// hgfPP.hpp
#ifndef HGFPP_H
#define HGFPP_H

#include <vector>
#include <string>

void
computeKTensorL ( const FluidMesh& Mesh, \
                  const std::vector<double>& xSolution, \
                  const std::vector<double>& ySolution, \
                  const std::vector<double>& zSolution, \
                  const std::vector<std::string>& KoutNames );

void
computeAveragesX ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution,
                   double& V, double& G, int print, \
                   const std::string& KoutName );

void
computeAveragesY ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print, \
                   const std::string& KoutName );

void
computeAveragesZ ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print, \
                   const std::string& KoutName );

void
computeKConstantDrive ( const FluidMesh& Mesh, \
                        const std::vector<double>& Solution,
                        double& K, \
                        int direction, \
                        int print, \
                        const std::string& KoutName );

void
computeKPoreNetwork( const PoreNetwork& pn, const std::vector<double>& Solution, \
                     const std::vector<double>& Ks, double& K, int direction, int print, \
                     const std::string& KoutName );

void
writeSolutionTP ( const FluidMesh& Mesh, const std::vector<double>& sol, \
                  std::string& outName );

void
writeSolutionPV ( const FluidMesh& Mesh, const std::vector<double>& sol, \
                  std::string& outName );

void
writePoreNetworkSolutionTP ( const PoreNetwork& pn, const std::vector<double>& sol, \
                             std::string& outName );

void
writeKCollection( const std::vector<double>& Kvec, int dim, std::string& outName );

#endif
