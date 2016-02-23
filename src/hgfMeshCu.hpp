// hgfMesh.hpp
#ifndef HGFMESH_H
#define HGFMESH_H

#include <vector>
#include <cuda_runtime.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "hgfAuxTools.hpp"

__global__ void ifcKernel2D ( unsigned long *d_CFC, const double *d_CCC, \
                              double epsx, double epsy, double xtol, double ytol, int nCells );
__global__ void ifcKernel3D ( unsigned long *d_CFC, const double *d_CCC, \
                              double epsx, double epsy, double epsz, \
                              double xtol, double ytol, double ztol, \
                              int nCells );
void MeshSubdivide( unsigned long *gridin, int ldi1, int ldi2, \
                    int nx, int ny, int nz, \
                    double length, double width, double height, \
                    int MX, int MY, int MZ, \
                    std::vector< std::vector<unsigned long> >& slices, \
                    std::vector<double>& lengths, \
                    std::vector<double>& widths, \
                    std::vector<double>& heights, \
                    std::vector<int>& nxs, \
                    std::vector<int>& nys, \
                    std::vector<int>& nzs );
void innerFaceConnectivity( std::vector<unsigned long>& ComponentFaceConnectivity, \
                            std::vector<double> ComponentCellCenters, \
                            double dx, double dy, double dz, int nCells, int DIM );

class FluidMesh
{
  public:
    // Public data
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
      ar & mv;
      ar & Nodes;
      ar & PCellCenters;
      ar & UCellCenters;
      ar & VCellCenters;
      ar & WCellCenters;
      ar & PCellWidths;
      ar & UCellWidths;
      ar & VCellWidths;
      ar & WCellWidths;
      ar & PresListByY;
      ar & PresListByZ;
      ar & PFaceConnectivity;
      ar & UFaceConnectivity;
      ar & VFaceConnectivity;
      ar & WFaceConnectivity;
      ar & PressureCellUNeighbor;
      ar & PressureCellVNeighbor;
      ar & PressureCellWNeighbor;
      ar & UCellPressureNeighbor;
      ar & VCellPressureNeighbor;
      ar & WCellPressureNeighbor;
      ar & UInteriorCells;
      ar & VInteriorCells;
      ar & WInteriorCells;
      ar & UBoundaryCells;
      ar & VBoundaryCells;
      ar & WBoundaryCells;
      ar & ImmersedBoundary;
      ar & FullGrid;
      ar & NodesLDI;
      ar & CellCentersLDI;
      ar & CellWidthsLDI;
      ar & FaceConnectivityLDI;
      ar & PressureCellVelocityNeighborLDI;
      ar & VelocityCellPressureNeighborLDI;
      ar & DIM;
      ar & NX;
      ar & NY;
      ar & NZ;
      ar & dofTotal;
      ar & maxNNZ;
      ar & xLim;
      ar & yLim;
      ar & zLim;
      ar & porosity;
    }
    std::vector<unsigned long> mv;
    std::vector<double> Nodes, PCellCenters, UCellCenters, VCellCenters;
    std::vector<double> WCellCenters;
    std::vector<double> PCellWidths, UCellWidths, VCellWidths, WCellWidths;
    std::vector<unsigned long> PresListByY, PresListByZ;
    std::vector<unsigned long> PFaceConnectivity, UFaceConnectivity;
    std::vector<unsigned long> VFaceConnectivity, WFaceConnectivity;
    std::vector<unsigned long> PressureCellUNeighbor, PressureCellVNeighbor;
    std::vector<unsigned long> PressureCellWNeighbor, UCellPressureNeighbor;
    std::vector<unsigned long> VCellPressureNeighbor, WCellPressureNeighbor;
    std::vector<unsigned long> UInteriorCells, VInteriorCells, WInteriorCells;
    std::vector<unsigned long> UBoundaryCells, VBoundaryCells, WBoundaryCells;
    std::vector<unsigned long> ImmersedBoundary, FullGrid;
    std::vector<int> DOF;
    int NodesLDI, CellCentersLDI, CellWidthsLDI, FaceConnectivityLDI;
    int PressureCellVelocityNeighborLDI, VelocityCellPressureNeighborLDI;
    int DIM, NX, NY, NZ, dofTotal, maxNNZ;
    double xLim[2], yLim[2], zLim[2], porosity;
    // Public functions
    int VelocityDOF( void );
    void BuildUniformMesh( const ProbParam& Par );
  private:
    int isNear2d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                  double dx, double dy, double dz, int nNodes );
    int isNear3d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                  double dx, double dy, double dz, int nNodes );
    void TotalDOF( void );
    void MaxNonZero( void );
    void sortPV( void );
};

class PoreNetwork
{
  public:
    // Public data
    std::vector<double> PoresXYZ;
    std::vector<unsigned long> Throats;
    std::vector<unsigned long> BoundaryPores, InteriorPores;
    int DIM, nPores, nThroats;
    double dx, dy, dz, psLength, psWidth, psHeight;
    // public functions
    void UniformPN( double length, double width, double height, int nx, int ny, int nz );
};

void SaveFluidMesh( const FluidMesh& Mesh, const std::string& outName );
void LoadFluidMesh( FluidMesh& Mesh, const std::string& inName );

#endif
