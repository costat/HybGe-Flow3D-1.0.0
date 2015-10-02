#include <vector>

void
ucartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction );

void
vcartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction );

void
wcartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction );

void
ucartflow2D ( const FluidMesh& Mesh, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals, \
              std::vector<double>& force, \
              double visc, int direction );

void
vcartflow2D ( const FluidMesh& Mesh, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals, \
              std::vector<double>& force, \
              double visc, int direction );

void
BoundaryConditions ( const FluidMesh& Mesh, double visc, \
                     std::vector<int>& matIs, std::vector<int>& matJs, \
                     std::vector<double>& matVals, \
                     std::vector<double>& force, int direction );
