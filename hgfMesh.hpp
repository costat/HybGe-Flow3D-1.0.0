// hgfMesh.hpp
#ifndef FluidMesh_H
#define FluidMesh_H

#include <vector>

class FluidMesh
{
  public:
    // Public data
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
    void BuildUniformMesh( unsigned long *gridin, int ldi1, int ldi2, \
                    int nx, int ny, int nz, \
                    double length, double width, double height );
  private:
    int isNear2d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                  double dx, double dy, double dz, int nNodes );
    int isNear3d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                  double dx, double dy, double dz, int nNodes );
    void innerFaceConnectivity( std::vector<unsigned long>& ComponentFaceConnectivity, \
                                std::vector<double> ComponentCellCenters, \
                                double dx, double dy, double dz, int nCells );
    void TotalDOF( void );
    void MaxNonZero( void );
    void sortPV( void );
/*  Some upcoming functions not yet implemented
    void AddCell( double *Nodes, int isIB );
    void RemoveCell( int GlobalCellNumber  );
    void IBGrowth( int GlobalCellNumber );
    void IBShrink( int GlobalCellNumber );
*/
};

#endif
