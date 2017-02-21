#!/usr/bin/env python

import numpy as np
import sys
import re
import os

if (len(sys.argv) != 5): # path not provided, quick exit
  print "Please provide geometry file, nx, ny, and nz for subdivision.\n"
  sys.exit(0)

dir_in = sys.argv[1]
nx_sd = int(sys.argv[2])
ny_sd = int(sys.argv[3])
nz_sd = int(sys.argv[4])

grid = open( dir_in+"/Geometry.dat", 'r' )
par = open( dir_in+"/Parameters.dat", 'r' )
grid.seek(0)
nx_data = grid.readline()
ny_data = grid.readline()
nz_data = grid.readline()
nx = re.findall(r'\d+', nx_data)
ny = re.findall(r'\d+', ny_data)
nz = re.findall(r'\d+', nz_data)
nx = int(nx[0])
ny = int(ny[0])
nz = int(nz[0])
grid_out = np.genfromtxt(grid, dtype = int) # import the rest as a numpy array
grid.close()
par.seek(0)
visc_data = par.readline()
visc = re.findall(r'\d+\.\d+', visc_data)
visc = float(visc[0])
length_data = par.readline()
length = re.findall(r'\d+\.\d+', length_data)
length = float(length[0])
width_data = par.readline()
width = re.findall(r'\d+\.\d+', width_data)
width = float(width[0])
height_data = par.readline()
height = re.findall(r'\d+\.\d+', height_data)
height = float(height[0])
par.close()

if nz:
  print "Creating ", nx_sd, "x", ny_sd, "x", nz_sd, " subdivision of ", nx, "x", ny, "x", nz, " mesh.\n"

  # set up subdivision sizes
  if (nx % nx_sd):
      sd_x_len = int(nx / nx_sd) + 1
  else:
      sd_x_len = int(nx / nx_sd)
  if (ny % ny_sd):
      sd_y_len = int(ny / ny_sd) + 1
  else:
      sd_y_len = int(ny / ny_sd)
  if (nz % nz_sd):
      sd_z_len = int(nz / nz_sd) + 1
  else:
      sd_z_len = int(nz / nz_sd)

  # loop through subdivisions
  for k in range(0,nz_sd):
      for j in range(0,ny_sd):
          for i in range(0,nx_sd):
              # create folder
              dir_name = '/Subdomain_'
              sd_num = i + (nx_sd * (j + ny_sd * k))
              dir_name += str(sd_num)
              if not os.path.exists(dir_in+dir_name):
                  os.makedirs(dir_in+dir_name)

              # create Parameters file
              pfilename = dir_in+dir_name+'/Parameters.dat'
              pfile = open(pfilename, 'w')
              pfile.write("viscosity= %f\n" % visc)
              sd_length = length * min(sd_x_len, nx-(i*sd_x_len)) / float(nx)
              pfile.write("length= %f\n" % sd_length)
              sd_width = width * min(sd_y_len, ny-(j*sd_y_len)) / float(ny)
              pfile.write("width= %f\n" % sd_width)
              sd_height = height * min(sd_z_len, nz-(k*sd_z_len)) / float(nz)
              pfile.write("height= %f\n" % sd_height)
              pfile.close()
              # create Geometry file
              gfilename = dir_in+dir_name+'/Geometry.dat'
              gfile = open(gfilename, 'w')
              gfile.write('nx= %d\n' % min(sd_x_len,nx-(i*sd_x_len)))
              gfile.write('ny= %d\n' % min(sd_y_len,ny-(j*sd_y_len)))
              gfile.write('nz= %d\n' % min(sd_z_len,nz-(k*sd_z_len)))
              for kk in range(k*sd_z_len,k*sd_z_len+min(sd_z_len,nz-(k*sd_z_len))):
                  for jj in range(j*sd_y_len,j*sd_y_len+min(sd_y_len,ny-(j*sd_y_len))):
                      for ii in range(i*sd_x_len,i*sd_x_len+min(sd_x_len,nx-(i*sd_x_len))):
                          gfile.write("%d " % grid_out[jj + ny * kk][ii])
                      gfile.write("\n")
                  gfile.write("\n")
              gfile.close()


else:
  print("Creating ", nx_sd, "x", ny_sd, " subdivision of ", nx, "x", ny, " mesh.\n" )

  # set up subdivision sizes
  if (nx % nx_sd):
      sd_x_len = int(nx / nx_sd) + 1
  else:
      sd_x_len = int(nx / nx_sd)
  if (ny % ny_sd):
      sd_y_len = int(ny / ny_sd) + 1
  else:
      sd_y_len = int(ny / ny_sd)

  # loop through subdivisions
  for j in range(0,ny_sd):
      for i in range(0,nx_sd):
          # create folder
          dir_name = '/Subdomain_'
          sd_num = i + nx_sd * j
          dir_name += str(sd_num)
          if not os.path.exists(dir_in+dir_name):
              os.makedirs(dir_in+dir_name)

          # create Parameters file
          pfilename = dir_in+dir_name+'/Parameters.dat'
          pfile = open(pfilename, 'w')
          pfile.write("viscosity= %f\n" % visc)
          sd_length = length * min(sd_x_len, nx-(i*sd_x_len)) / float(nx)
          pfile.write("length= %f\n" % sd_length)
          sd_width = width * min(sd_y_len, ny-(j*sd_y_len)) / float(ny)
          pfile.write("width= %f\n" % sd_width)
          sd_height = 0
          pfile.write("height= %f\n" % sd_height)
          pfile.close()
          # create Geometry file
          gfilename = dir_in+dir_name+'/Geometry.dat'
          gfile = open(gfilename, 'w')
          gfile.write('nx= %d\n' % min(sd_x_len,nx-(i*sd_x_len)))
          gfile.write('ny= %d\n' % min(sd_y_len,ny-(j*sd_y_len)))
          gfile.write('nz= %d\n' % 0)
          for jj in range(j*sd_y_len,j*sd_y_len+min(sd_y_len,ny-(j*sd_y_len))):
              for ii in range(i*sd_x_len,i*sd_x_len+min(sd_x_len,nx-(i*sd_x_len))):
                  gfile.write("%d " % grid_out[jj][ii])
              gfile.write("\n")
          gfile.close()
