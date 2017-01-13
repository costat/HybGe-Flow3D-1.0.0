/* utilities source */

// system includes
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// 1d->3d index
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

#include "../../include/hgflow.hpp"

namespace bfs = boost::filesystem;

void
hgf::init_parameters(parameters& par, const std::string& problem_path)
{
  par.problem_path = problem_path;
  std::string param = "Parameters.dat";
  bfs::path Parameters;
  bool isParam = hgf::find_file(par.problem_path, param, Parameters);
  if (!isParam) {
	  std::cout << "\nParameter file not present in problem folder. Exiting\n";
	  exit(0);
  }
  hgf::load_parameters(par, Parameters);
  hgf::import_voxel_geometry(par, par.problem_path);
  par.verbose = 0;
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
  std::cout << "Dimension= " << par.dimension << "\n";
  std::cout << "Geometry length= " << par.length << "\n";
  std::cout << "Geometry width= " << par.width << "\n";
  std::cout << "Geometry height= " << par.height << "\n";
  std::cout << "Verbose= " << par.verbose << "\n";
  std::cout << "Problem path= " << par.problem_path.string() << "\n";
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

// function removes cells that are boundaries in opposite directions
int 
hgf::mesh_sanity(parameters& par)
{
  int totalChanged = 0;
  int nChanged;

  if (par.dimension == 3) goto sanityCheck3;
  else goto sanityCheck2;

sanityCheck3:
  {
    nChanged = 0;
    // sanity
    for (int zi = 0; zi < par.nz; zi++) {
      for (int yi = 0; yi < par.ny; yi++) {
        for (int xi = 0; xi < par.nx; xi++) {
          if (par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] != 1) {
            // xi sanity
            if (xi == 0) {
              if (par.voxel_geometry[idx3(zi, yi, (xi + 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (xi == par.nx - 1) {
              if (par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3(zi, yi, (xi + 1), par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            // yi sanity
            if (yi == 0) {
              if (par.voxel_geometry[idx3(zi, (yi + 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (yi == par.ny - 1) {
              if (par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3(zi, (yi + 1), xi, par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            // zi sanity
            if (zi == 0) {
              if (par.voxel_geometry[idx3((zi + 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (zi == par.nz - 1) {
              if (par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3((zi + 1), yi, xi, par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck3;
    else goto cleanup;
  }

sanityCheck2:
  {
    nChanged = 0;
    // sanity
    for (int yi = 0; yi < par.ny; yi++) {
      for (int xi = 0; xi < par.nx; xi++) {
        if (par.voxel_geometry[idx2(yi, xi, par.nx)] != 1) {
          // xi sanity
          if (xi == 0) {
            if (par.voxel_geometry[idx2(yi, (xi + 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else if (xi == par.nx - 1) {
            if (par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else {
            if (par.voxel_geometry[idx2(yi, (xi + 1), par.nx)] == 1 \
              && par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          // yi sanity
          if (yi == 0) {
            if (par.voxel_geometry[idx2((yi + 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else if (yi == par.ny - 1) {
            if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else {
            if (par.voxel_geometry[idx2((yi + 1), xi, par.nx)] == 1 \
              && par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck2;
    else goto cleanup;
  }

cleanup:
  if (totalChanged) {
    std::cout << "\nWarning, input geometry was incompatible.\n";
    std::cout << totalChanged << " cells, representing ";
    if (par.dimension == 3) std::cout << (double)100 * totalChanged / (par.nx*par.ny*par.nz);
    else std::cout << (double)100 * totalChanged / (par.nx*par.ny);
    std::cout << "% of the input geometry, with boundaries on opposite faces \nwere found and removed from void space.\n";
  }
  return totalChanged;
}

int
hgf::remove_dead_pores(parameters& par)
{
  std::vector<unsigned long> voxel_geometry_cpy(par.voxel_geometry);
  std::vector<unsigned long> search_queue;
  std::vector<int> current_component;
  int n_components = 0;
  int pores_removed = 0;
  int ii, jj, kk;
  if (par.dimension == 3) {
    for (int k = 0; k < par.nz; k++) {
      for (int j = 0; j < par.ny; j++) {
        for (int i = 0; i < par.nx; i++) {
          if (voxel_geometry_cpy[idx3(k, j, i, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(k, j, i, par.ny, par.nx)] == 2) {
            ii = i;
            jj = j;
            kk = k;
            n_components++;
            voxel_geometry_cpy[idx3(kk, jj, ii, par.ny, par.nx)] = n_components + 2;
            current_component.resize(3);
            current_component[0] = k;
            current_component[1] = j;
            current_component[2] = i;
          }
          else continue;
          search_3d:
          // look y-
          if (jj) {
            if (voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back((jj - 1));
              current_component.push_back(ii);
              search_queue.push_back(kk);
              search_queue.push_back((jj - 1));
              search_queue.push_back(ii);
            }
          }
          // look x+
          if (ii < par.nx - 1) {
            if (voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back(jj);
              current_component.push_back((ii + 1));
              search_queue.push_back(kk);
              search_queue.push_back(jj);
              search_queue.push_back((ii + 1));
            }
          }
          // look y+
          if (jj < par.ny - 1) {
            if (voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back((jj + 1));
              current_component.push_back(ii);
              search_queue.push_back(kk);
              search_queue.push_back((jj + 1));
              search_queue.push_back(ii);
            }
          }
          // look x-
          if (ii) {
            if (voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back(jj);
              current_component.push_back((ii - 1));
              search_queue.push_back(kk);
              search_queue.push_back(jj);
              search_queue.push_back((ii - 1));
            }
          }
         
          // look z-
          if (kk) {
            if (voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back((kk - 1));
              current_component.push_back(jj);
              current_component.push_back(ii);
              search_queue.push_back((kk - 1));
              search_queue.push_back(jj);
              search_queue.push_back(ii);
            }
          }
          // look z+
          if (kk < par.nz - 1) {
            if (voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back((kk + 1));
              current_component.push_back(jj);
              current_component.push_back(ii);
              search_queue.push_back((kk + 1));
              search_queue.push_back(jj);
              search_queue.push_back(ii);
            }
          }
          if (search_queue.size()) {
            ii = search_queue[search_queue.size() - 1];
            jj = search_queue[search_queue.size() - 2];
            kk = search_queue[search_queue.size() - 3];
            search_queue.resize(search_queue.size() - 3);
            goto search_3d;
          }
          else { // done filling in this component, check if pore needs removing
            int min_i = par.nx;
            int min_j = par.ny;
            int min_k = par.nz;
            int max_i = 0;
            int max_j = 0;
            int max_k = 0;
            for (int cc = 0; cc < ((int)current_component.size() / 3); cc++) {
              if (current_component[idx2(cc, 0, 3)] > max_k) max_k = current_component[idx2(cc, 0, 3)];
              if (current_component[idx2(cc, 0, 3)] < min_k) min_k = current_component[idx2(cc, 0, 3)];
              if (current_component[idx2(cc, 1, 3)] > max_j) max_j = current_component[idx2(cc, 1, 3)];
              if (current_component[idx2(cc, 1, 3)] < min_j) min_j = current_component[idx2(cc, 1, 3)];
              if (current_component[idx2(cc, 2, 3)] > max_i) max_i = current_component[idx2(cc, 2, 3)];
              if (current_component[idx2(cc, 2, 3)] < min_i) min_i = current_component[idx2(cc, 2, 3)];
            }
            if (min_i > 0 || min_j > 0 || max_i < par.nx - 1 || max_j < par.ny - 1 || min_k > 0 || max_k < par.nz - 1) {
              pores_removed++;
              for (int cc = 0; cc < ((int)current_component.size() / 3); cc++) {
                par.voxel_geometry[idx3(current_component[idx2(cc, 0, 3)], current_component[idx2(cc, 1, 3)], current_component[idx2(cc, 2, 3)], par.ny, par.nx)] = 1;
              }
            }
          }
        }
      }
    }
  }
  else {
    for (int j = 0; j < par.ny; j++) {
      for (int i = 0; i < par.nx; i++) {
        if (voxel_geometry_cpy[idx2(j, i, par.nx)] == 0 || voxel_geometry_cpy[idx2(j, i, par.nx)] == 2) {
          ii = i;
          jj = j;
          n_components++;
          voxel_geometry_cpy[idx2(jj, ii, par.nx)] = n_components + 2;
          current_component.resize(2);
          current_component[0] = j;
          current_component[1] = i;
        }
        else continue;
        search_2d:
        // look down
        if (jj) {
          if (voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] == 0 || voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] == 2) {
            voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] = n_components + 2;
            current_component.push_back((jj - 1));
            current_component.push_back(ii);
            search_queue.push_back((jj - 1));
            search_queue.push_back(ii);
          }
        }
        // look right
        if (ii < par.nx - 1) {
          if (voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] == 0 || voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] == 2) {
            voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] = n_components + 2;
            current_component.push_back(jj);
            current_component.push_back((ii + 1));
            search_queue.push_back(jj);
            search_queue.push_back((ii + 1));
          }
        }
        // look up
        if (jj < par.ny - 1) {
          if (voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] == 0 || voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] == 2) {
            voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] = n_components + 2;
            current_component.push_back((jj + 1));
            current_component.push_back(ii);
            search_queue.push_back((jj + 1));
            search_queue.push_back(ii);
          }
        }
        // look left 
        if (ii) {
          if (voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] == 0 || voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] == 2) {
            voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] = n_components + 2;
            current_component.push_back(jj);
            current_component.push_back((ii - 1));
            search_queue.push_back(jj);
            search_queue.push_back((ii - 1));
          }
        }
        if (search_queue.size()) {
          ii = search_queue[search_queue.size() - 1];
          jj = search_queue[search_queue.size() - 2];
          search_queue.resize(search_queue.size() - 2);
          goto search_2d;
        }
        else { // done filling in this component, check if pore needs removing
          int min_i = par.nx;
          int min_j = par.ny;
          int max_i = 0;
          int max_j = 0;
          for (int cc = 0; cc < ((int)current_component.size() / 2); cc++) {
            if (current_component[idx2(cc, 0, 2)] > max_j) max_j = current_component[idx2(cc, 0, 2)];
            if (current_component[idx2(cc, 0, 2)] < min_j) min_j = current_component[idx2(cc, 0, 2)];
            if (current_component[idx2(cc, 1, 2)] > max_i) max_i = current_component[idx2(cc, 1, 2)];
            if (current_component[idx2(cc, 1, 2)] < min_i) min_i = current_component[idx2(cc, 1, 2)];
          }
          if (min_i > 0 || min_j > 0 || max_i < par.nx - 1 || max_j < par.ny - 1) {
            pores_removed++;
            for (int cc = 0; cc < ((int)current_component.size() / 2); cc++) {
              par.voxel_geometry[idx2(current_component[idx2(cc, 0, 2)], current_component[idx2(cc, 1, 2)], par.nx)] = 1;
            }
          }
        }
      }
    }
  }

  return pores_removed;
}
