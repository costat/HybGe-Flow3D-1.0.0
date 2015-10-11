## Script to compute constant Ks (x direction) for all geometries in a directory ###

import numpy as np
import sys
import hgf
import re

gridfolder = '/media/basis1/global/research/hybrid/kfactory/ranThroats/'
L = 1.
W = 1.
H = 1.

direction = 0
visc = 0.001
nThreads = 4

## Cycle through grids in directory ##
for root, dirs, filenames in os.walk(gridfolder):
  for f in filenames:
    gridF = open(os.path.join(root, f), 'r')
    gridF.seek(0)
    geoData = gridF.readline()
    gridin1 = np.genfromtxt(gridF, dtype = int)
    gridF.close()
    nxyz = re.findall(r'\d+', geoData)
    nx = int(nxzy[0])
    ny = int(nxyz[1])
    nz = int(nxyz[2])

    if nz:
      gridin = np.zeros((nx, ny, nz), dtype = int)
      for zLevel in range(0,nz):
        gridin[:,:,zLevel] = gridin1[zLevel*ny:(zLevel+1)*ny][:]
      gridin_ldi1, gridin_ldi2, gridin_ldi3 = gridin.shape
      gridin = np.array(gridin, dtype = np.uint64).ravel()
    elif not nz:
      gridin = np.transpose(gridin1)
      gridin_ldi1, gridin_ldi2 = gridin.shape
      gridin = np.array(gridin, dtype = np.uint64).ravel()
      gridin_ldi3 = 0

    gridin_len = gridin.size

    print 'Solving the stationary Stokes problem...\n'

    #################
    ### SWIG CALL ###
    #################

    hgf.hgfStokesDrive ( gridin, gridin_ldi2, gridin_ldi3, nx, ny, nz, \
                     L, W, H, direction, visc, nThreads )

    ### Grab computed K and append to file to hold all Ks ###

    ### Move flow solution file so not overwritten ###

