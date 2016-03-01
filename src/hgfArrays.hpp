// hgfArrays.hpp
#ifndef HGFARRAYS_H
#define HGFARRAYS_H

#include <vector>
#include "hgfAuxTools.hpp"

void
diffusionArrays ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  const std::vector<unsigned long>& velint, \
                  const std::vector<unsigned long>& fconn, \
                  const std::vector<double>& dxyz, \
                  double visc, int velShift );

void
diffusionDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, \
                  std::vector<double>& matVals, \
                  int Shift, double visc );

void
pressureGrad ( const FluidMesh& Mesh, std::vector<int>& matIs, \
               std::vector<int>& matJs, std::vector<double>& matVals, \
               const std::vector<unsigned long>& velint, \
               const std::vector<unsigned long>& vcpn, \
               const std::vector<double>& dxyz, \
               int velShift, int done, int dtwo );

void
pressureGradDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                    std::vector<int>& matJs, \
                    std::vector<double>& matVals, \
                    int Shift );

void
continuityEq ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                std::vector<int>& matJs, std::vector<double>& matVals, \
                const std::vector<unsigned long>& pcvn, \
                int velShift, int done, int dtwo );

void
continuityDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, \
                  std::vector<double>& matVals, \
                  int Shift );

void
GradientTranspose( const FluidMesh& Mesh, std::vector<int>& matIs, \
                   std::vector<int>& matJs, std::vector<double>& matVals );

void
StokesArray ( const FluidMesh& Mesh, double visc, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals );

void
VelocityArray ( const FluidMesh& Mesh, double visc, std::vector<int>& matIs, \
                std::vector<int>& matJs, std::vector<double>& matVals, int component );

void
PoreNetworkArray ( const PoreNetwork& pn, std::vector<int>& matIs, \
                   std::vector<int>& matJs, std::vector<double>& matVals, \
                   const std::vector<double>& Ks );

void
sortCOO ( std::vector<int>& matIs, std::vector<int>& matJs, std::vector<double>& matVals );

void
buildCSR ( const std::vector<int>& matIs, const std::vector<int>& matJs, \
           const std::vector<double>& matVals, std::vector<int>& rowPTR );

#endif
