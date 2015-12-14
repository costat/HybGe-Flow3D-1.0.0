#include <vector>

void
FlowComponent3d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, \
                  const std::vector<unsigned long>& componentBoundary, \
                  const std::vector<double>& componentCellCenters, \
                  const std::vector<double>& componentCellWidths, \
                  const std::vector<double>& CellWidthsLR, \
                  const std::vector<double>& CellWidthsUD, \
                  const std::vector<unsigned long>& componentConnectivity, \
                  const std::vector<unsigned long>& PressureNeighbor, \
                  double visc, int direction, int component );

void
FlowComponent2d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, \
                  const std::vector<unsigned long>& componentBoundary, \
                  const std::vector<double>& componentCellCenters, \
                  const std::vector<double>& componentCellWidths, \
                  const std::vector<double>& CellWidthsLR, \
                  const std::vector<unsigned long>& componentConnectivity, \
                  const std::vector<unsigned long>& PressureNeighbor, \
                  double visc, int direction, int component );

void
AxisFlowDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                std::vector<int>& matJs, std::vector<double>& matVals, \
                std::vector<double>& force, double visc, int direction );
