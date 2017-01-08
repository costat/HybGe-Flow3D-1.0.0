#include "multiscale_flow.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

double
hgf::multiscale::flow::compute_permeability_x(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                                     const std::vector< degree_of_freedom >& velocity_w, \
                                                                     const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  double min_x, max_x, mid_x, midrange_x, min_y, max_y, min_z, max_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  if (par.dimension == 3) {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    min_z = 0;
    max_z = par.height;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    mid_x = 0.5 * (min_x + max_x);
    min_y += r*par.width;
    max_y -= r*par.width;
    min_z += r*par.height;
    max_z -= r*par.height;
    midrange_x = 0.5 * (max_x + mid_x) - 0.5 * (mid_x + min_x);
    // compute averages
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      p1_idx = velocity_u[ii].cell_numbers[0];
      p2_idx = velocity_u[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1) {
        if (velocity_u[ii].coords[0] > min_x) {
          if (velocity_u[ii].coords[0] < max_x) {
            if (velocity_u[ii].coords[1] > min_y) {
              if (velocity_u[ii].coords[1] < max_y) {
                if (velocity_u[ii].coords[2] > min_z) {
                  if (velocity_u[ii].coords[2] < max_z) {
                    pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                    if (velocity_u[ii].coords[0] < mid_x) {
                      p1 += pressure;
                      p1_count++;
                      v += solution[ii];
                      v_count++;
                    }
                    else if (velocity_u[ii].coords[0] >= mid_x) {
                      p2 += pressure;
                      p2_count++;
                      v += solution[ii];
                      v_count++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // compute v, g 2d
  else {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    mid_x = 0.5 * (min_x + max_x);
    min_y += r*par.width;
    max_y -= r*par.width;
    midrange_x = 0.5 * (max_x + mid_x) - 0.5 * (mid_x + min_x);
    // compute averages
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      p1_idx = velocity_u[ii].cell_numbers[0];
      p2_idx = velocity_u[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1) {
        if (velocity_u[ii].coords[0] > min_x) {
          if (velocity_u[ii].coords[0] < max_x) {
            if (velocity_u[ii].coords[1] > min_y) {
              if (velocity_u[ii].coords[1] < max_y) {
                pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                if (velocity_u[ii].coords[0] < mid_x) {
                  p1 += pressure;
                  p1_count++;
                  v += solution[ii];
                  v_count++;
                }
                else if (velocity_u[ii].coords[0] >= mid_x) {
                  p2 += pressure;
                  p2_count++;
                  v += solution[ii];
                  v_count++;
                }
              }
            }
          }
        }
      }
    }
  }
  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_x;

  // determine porosity (holder)
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v/g) * por;

}

double 
hgf::multiscale::flow::compute_permeability_y(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                     const std::vector< degree_of_freedom >& velocity_w, \
  const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  double min_x, max_x, min_y, max_y, mid_y, midrange_y, min_z, max_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  if (par.dimension == 3) {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    min_z = 0;
    max_z = par.height;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    min_y += r*par.width;
    max_y -= r*par.width;
    mid_y = 0.5 * (min_y + max_y);
    min_z += r*par.height;
    max_z -= r*par.height;
    midrange_y = 0.5 * (max_y + mid_y) - 0.5 * (mid_y + min_y);
    // compute averages
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      p1_idx = velocity_v[ii].cell_numbers[0];
      p2_idx = velocity_v[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1) {
        if (velocity_v[ii].coords[0] > min_x) {
          if (velocity_v[ii].coords[0] < max_x) {
            if (velocity_v[ii].coords[1] > min_y) {
              if (velocity_v[ii].coords[1] < max_y) {
                if (velocity_v[ii].coords[2] > min_z) {
                  if (velocity_v[ii].coords[2] < max_z) {
                    pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                    if (velocity_v[ii].coords[0] < mid_y) {
                      p1 += pressure;
                      p1_count++;
                      v += solution[velocity_u.size() + ii];
                      v_count++;
                    }
                    else if (velocity_v[ii].coords[0] >= mid_y) {
                      p2 += pressure;
                      p2_count++;
                      v += solution[velocity_u.size() + ii];
                      v_count++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // compute v, g 2d
  else {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    min_y += r*par.width;
    max_y -= r*par.width;
    mid_y = 0.5 * (min_y + max_y);
    midrange_y = 0.5 * (max_y + mid_y) - 0.5 * (mid_y + min_y);
    // compute averages
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      p1_idx = velocity_v[ii].cell_numbers[0];
      p2_idx = velocity_v[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1) {
        if (velocity_v[ii].coords[0] > min_x) {
          if (velocity_v[ii].coords[0] < max_x) {
            if (velocity_v[ii].coords[1] > min_y) {
              if (velocity_v[ii].coords[1] < max_y) {
                pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                if (velocity_v[ii].coords[0] < mid_y) {
                  p1 += pressure;
                  p1_count++;
                  v += solution[velocity_u.size() + ii];
                  v_count++;
                }
                else if (velocity_v[ii].coords[0] >= mid_y) {
                  p2 += pressure;
                  p2_count++;
                  v += solution[velocity_u.size() + ii];
                  v_count++;
                }
              }
            }
          }
        }
      }
    }
  }
  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_y;

  // determine porosity (holder)
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v / g) * por;

}

double 
hgf::multiscale::flow::compute_permeability_z(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                     const std::vector< degree_of_freedom >& velocity_w, \
                                                     const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  double min_x, max_x, min_y, max_y, min_z, max_z, mid_z, midrange_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
  // geometry properties
  min_x = 0;
  max_x = par.length;
  min_y = 0;
  max_y = par.width;
  min_z = 0;
  max_z = par.height;
  // adjust limits to avoid recirculation
  min_x += r*par.length;
  max_x -= r*par.length;
  min_y += r*par.width;
  max_y -= r*par.width;
  min_z += r*par.height;
  max_z -= r*par.height;
  mid_z = 0.5 * (min_z + max_z);
  midrange_z = 0.5 * (max_z + mid_z) - 0.5 * (mid_z + min_z);
  // compute averages
  for (int ii = 0; ii < velocity_w.size(); ii++) {
    p1_idx = velocity_w[ii].cell_numbers[0];
    p2_idx = velocity_w[ii].cell_numbers[1];
    if (p1_idx != -1 && p2_idx != -1) {
      if (velocity_w[ii].coords[0] > min_x) {
        if (velocity_w[ii].coords[0] < max_x) {
          if (velocity_w[ii].coords[1] > min_y) {
            if (velocity_w[ii].coords[1] < max_y) {
              if (velocity_w[ii].coords[2] > min_z) {
                if (velocity_w[ii].coords[2] < max_z) {
                  pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                  if (velocity_w[ii].coords[0] < mid_z) {
                    p1 += pressure;
                    p1_count++;
                    v += solution[velocity_u.size() + velocity_v.size() + ii];
                    v_count++;
                  }
                  else if (velocity_w[ii].coords[0] >= mid_z) {
                    p2 += pressure;
                    p2_count++;
                    v += solution[velocity_u.size() + velocity_v.size() + ii];
                    v_count++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_z;

  // determine porosity (holder)
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v / g) * por;
}
