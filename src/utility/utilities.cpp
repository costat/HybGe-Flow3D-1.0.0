/* utilities source */

// system includes
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>

#include "../../include/hgflow.hpp"

namespace bfs = boost::filesystem;

void
hgf::init_parameters(parameters& par, const std::string& problem_path)
{
  par.problem_path = problem_path;
  std::cout << "\nSolving problem in directory: \t " << par.problem_path << "\n";
  std::string param = "Parameters.dat";
  bfs::path Parameters;
  bool isParam = hgf::find_file(par.problem_path, param, Parameters);
  if (!isParam) {
	  std::cout << "\nParameter file not present in problem folder. Exiting\n";
	  exit(0);
  }
  hgf::load_parameters(par, Parameters);
  hgf::import_voxel_geometry(par, par.problem_path);
  hgf::print_parameters(par);
}

bool
hgf::find_file(const bfs::path& problem_path, \
  const std::string& file_name, \
  bfs::path& file_path)
{
  if (!exists(problem_path)) return false;
  bfs::directory_iterator end_itr;
  for (bfs::directory_iterator itr(problem_path); itr != end_itr; ++itr) {
    if (bfs::is_directory(itr->status())) {
      if (hgf::find_file(itr->path(), file_name, file_path)) return true;
    }
    else if (itr->path().leaf() == file_name) {
      file_path = itr->path();
      return true;
    }
  }
  return false;
}

void
hgf::load_parameters(parameters& par, const bfs::path& problem_path)
{
  std::string line;
  std::string str;
  bfs::ifstream ifs(problem_path);

  //--- viscosity ---//
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iVisc(line);
  iVisc >> str >> par.viscosity;

  //--- bc ---//
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iBC(line);
  iBC >> str >> par.bc;

  //--- length ---//
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ilength(line);
  ilength >> str >> par.length;

  //--- width ---//
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iwidth(line);
  iwidth >> str >> par.width;

  //--- height ---//
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iheight(line);
  iheight >> str >> par.height;
}

void
hgf::print_parameters(parameters& par)
{
  std::cout << "Viscosity= " << par.viscosity << "\n";
  std::cout << "Boundary Condtions= " << par.bc << "\n";
  std::cout << "Dimension= " << par.dimension << "\n";
  std::cout << "Geometry length= " << par.length << "\n";
  std::cout << "Geometry width= " << par.width << "\n";
  std::cout << "Geometry height= " << par.height << "\n";
}

void
hgf::import_voxel_geometry(parameters& par, const bfs::path& problem_path)
{
  std::string geometry_file = "Geometry.dat";
  bfs::path geo;
  bool isGeo = find_file(problem_path, geometry_file, geo);

  // error and exit if geometry file is missing
  if (!isGeo) {
    std::cout << "\nGeometry file not present and voxel geometry was requested. Exiting.\n";
    exit(0);
  }

  std::string line;
  std::string str;
  bfs::ifstream ifs(geo);

  // grab nx
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ix(line);
  ix >> str >> par.nx;

  // grab ny
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iy(line);
  iy >> str >> par.ny;

  // grab nz
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iz(line);
  iz >> str >> par.nz;

  if (!par.nz) par.dimension = 2;
  else par.dimension = 3;

  if (!par.nx || !par.ny) {
    std::cout << "\nSomething went wrong reading mesh. Ensure correct formatting. Exiting.\n";
    exit(0);
  }

  std::string hold;
  // read remaining lines of geometry into parameters file
  if (par.nz) { // 3d voxel file
    for (int nslices = 0; nslices < par.nz; nslices++) {
      for (int nrows = 0; nrows < par.ny; nrows++) {
        if (ifs.good()) {
          std::getline(ifs, line);
        }
        // place line entries into parameters data
        std::stringstream stream(line);
        while (1) {
          int n;
          stream >> n;
          if (!stream) break;
          par.voxel_geometry.push_back(n);
        }
      }
      // add the empty row between z slices
      if (ifs.good()) {
        std::getline(ifs, line);
      }
    }
  }
  else { // 2d voxel file
    for (int nrows = 0; nrows < par.ny; nrows++) {
      // grab row
      if (ifs.good()) {
        std::getline(ifs, line);
      }
      // place line entries into param data
      std::stringstream stream(line);
      while (1) {
        int n;
        stream >> n;
        if (!stream) break;
        par.voxel_geometry.push_back(n);
      }
    }
  }

}

bool
hgf::check_symmetry(std::vector< array_coo >& array)
{

  std::vector< array_coo > temp_array(array);
  hgf::unique_array(temp_array);
  std::vector< array_coo > array_trans(temp_array);
  std::sort(array_trans.begin(), array_trans.end(), byJbyI());

  for (int ii = 0; ii < temp_array.size(); ii++) {
    if ((temp_array[ii].i_index != array_trans[ii].j_index) || (temp_array[ii].j_index != array_trans[ii].i_index) || fabs(temp_array[ii].value - array_trans[ii].value) > 1e-12) {
      return 0;
    }
  }
  return 1;

}

// sorts coo array by rows then by columns (j runs min to max within 1 i)
void
hgf::sort_array(std::vector< array_coo >& array)
{

  std::sort(array.begin(), array.end(), byIbyJ());

 }

// create unique coo array. array sorted and then values associated to repeat indexes are summed, and extras removed
void
hgf::unique_array(std::vector< array_coo >& array)
{

  std::sort(array.begin(), array.end(), byIbyJ());

  int nthrs = omp_get_max_threads();
  int block_size = ((int)array.size() % nthrs) ? ((int)array.size() / nthrs) : ((int)array.size() / nthrs + 1);
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(nthrs);
  for (int ii = 0; ii < nthrs; ii++) temp_arrays[ii].reserve(block_size);

#pragma omp parallel shared(array)
  {
#pragma omp for schedule(static,1)
    for (int ii = 0; ii < nthrs; ii++) {
      int cpy_beg = ii * block_size;
      int cpy_end = cpy_beg;
      int cpy_size = 0;
      int icount = cpy_beg + 1;
      int tar_beg = 0;
      do {
        if (array[icount].i_index == array[icount - 1].i_index && array[icount].j_index == array[icount - 1].j_index) {
          cpy_end = icount - 1;
          array[cpy_end].value += array[icount].value;
          icount++;
          while (array[icount].i_index == array[cpy_end].i_index && array[icount].j_index == array[cpy_end].j_index && icount < std::min((ii+1)*block_size, (int)array.size())) {
            array[cpy_end].value += array[icount].value;
            icount++;
          }
          // copy block to temp
          cpy_size = cpy_end - cpy_beg + 1;
          temp_arrays[ii].resize(temp_arrays[ii].size() + cpy_size);
          memcpy(temp_arrays[ii].data() + tar_beg, array.data() + cpy_beg, cpy_size * sizeof(array_coo));

          // set idx for next iter
          cpy_beg = icount;
          tar_beg += cpy_size;
          icount++;
        }
        else {
          icount++;
        }
      } while (icount < std::min((ii+1)*block_size, (int)array.size()));
      if (cpy_beg < std::min((ii + 1)*block_size, (int)array.size())) {
        // final copy
        cpy_size = std::min((ii + 1)*block_size, (int)array.size()) - cpy_beg;
        temp_arrays[ii].resize(temp_arrays[ii].size() + cpy_size);
        memcpy(temp_arrays[ii].data() + tar_beg, array.data() + cpy_beg, cpy_size * sizeof(array_coo));
      }
    }
  }
  // paste
  int cpy_idx;
  memcpy(array.data(), temp_arrays[0].data(), temp_arrays[0].size() * sizeof(array_coo));
  cpy_idx = (int)temp_arrays[0].size();
  for (int ii = 1; ii < nthrs; ii++) {
    if (temp_arrays[ii - 1][temp_arrays[ii - 1].size() - 1].i_index == temp_arrays[ii][0].i_index && temp_arrays[ii - 1][temp_arrays[ii - 1].size() - 1].j_index == temp_arrays[ii][0].j_index) {
      cpy_idx--;
    }
    memcpy(array.data() + cpy_idx, temp_arrays[ii].data(), temp_arrays[ii].size() * sizeof(array_coo));
    cpy_idx += (int)temp_arrays[ii].size();
  }
  array.resize(cpy_idx);
  
}
