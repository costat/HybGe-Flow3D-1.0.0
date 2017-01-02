/* utilities header */
#ifndef _UTILITIES_H
#define _UTILITIES_H

// system includes
#include <vector>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

namespace hgf
{
  void
  init_parameters(parameters& par, const std::string& problem_path);
  
  bool
  find_file( const bfs::path& problem_path, \
             const std::string& file_name, \
             bfs::path& file_path);

  void
  load_parameters(parameters& par, const bfs::path& problem_path);

  void
  print_parameters(parameters& par);

  void
  import_voxel_geometry(parameters& par, const bfs::path& problem_path);

  bool
  check_symmetry(array_coo& arr);

  void
  unique_array(array_coo& arr);
}

#endif
