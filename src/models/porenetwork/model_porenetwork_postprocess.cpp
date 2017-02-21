/* porenetwork postprocess source */

// hgf includes
#include "model_porenetwork.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::uniform_porenetwork::output_vtk saves the solution to the porenetwork flow to a file for VTK visualiztion.
 *
 * @param[in] par - parameters struct containing problem information, including problem directory.
 */
void
hgf::models::uniform_porenetwork::output_vtk(const parameters& par)
{
  int nlines = 0;
  for (int ii = 0; ii < (int)pressure.size(); ii++) {
    for (int jj = 0; jj < (par.dimension * 2); jj++) {
      if (pressure[ii].neighbors[jj] != -1) nlines++;
    }
  }
  nlines /= 2;


  if (par.dimension == 3) {

    // write to vtk file
    bfs::path output_path(par.problem_path / "PN_Solution.vtk");
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET POLYDATA\n";
    outstream << "POINTS " << pressure.size() << " double\n";

    double lvalue;
    for (int row = 0; row < (int)pressure.size(); row++) {
      outstream << pressure[row].coords[0] << "\t";
      outstream << pressure[row].coords[1] << "\t";
      outstream << pressure[row].coords[2] << "\n";
    }
    outstream << "\n";
    outstream << "LINES " << nlines << " " << nlines * 3 << "\n";
    for (int row = 0; row < (int)pressure.size(); row++) {
      for (int nbr = 0; nbr < 6; nbr++) {
        if (pressure[row].neighbors[nbr] > row) {
          outstream << 2 << "\t" << row << "\t" << pressure[row].neighbors[nbr] << "\n";
        }
      }
    }
    outstream << "\n";
    outstream << "POINT_DATA " << pressure.size() << "\n";
    outstream << "SCALARS pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int ii = 0; ii < (int)pressure.size(); ii++) {
      outstream << solution[ii] << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nlines << "\n";
    outstream << "SCALARS throat_pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < (int)pressure.size(); row++) {
      for (int nbr = 0; nbr < 6; nbr++) {
        if (pressure[row].neighbors[nbr] > row) {
          lvalue = 0.5*(solution[row] + solution[pressure[row].neighbors[nbr]]);
          outstream << lvalue << "\n";
        }
      }
    }
    outstream << "\n";

    outstream.close();

  }
  else {

    // write to vtk file
    bfs::path output_path(par.problem_path / "PN_Solution.vtk");
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET POLYDATA\n";
    outstream << "POINTS " << pressure.size() << " double\n";

    double lvalue;
    for (int row = 0; row < (int)pressure.size(); row++) {
      outstream << pressure[row].coords[0] << "\t";
      outstream << pressure[row].coords[1] << "\t";
      outstream << 0.0 << "\n";
    }
    outstream << "\n";
    outstream << "LINES " << nlines << " " << nlines * 3 << "\n";
    for (int row = 0; row < (int)pressure.size(); row++) {
      for (int nbr = 0; nbr < 4; nbr++) {
        if (pressure[row].neighbors[nbr] > row) {
          outstream << 2 << "\t" << row << "\t" << pressure[row].neighbors[nbr] << "\n";
        }
      }
    }
    outstream << "\n";
    outstream << "POINT_DATA " << pressure.size() << "\n";
    outstream << "SCALARS pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int ii = 0; ii < (int)pressure.size(); ii++) {
      outstream << solution[ii] << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nlines << "\n";
    outstream << "SCALARS throat_pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < (int)pressure.size(); row++) {
      for (int nbr = 0; nbr < 4; nbr++) {
        if (pressure[row].neighbors[nbr] > row) {
          lvalue = 0.5*(solution[row] + solution[pressure[row].neighbors[nbr]]);
          outstream << lvalue << "\n";
        }
      }
    }
    outstream << "\n";

    outstream.close();

  }

}
