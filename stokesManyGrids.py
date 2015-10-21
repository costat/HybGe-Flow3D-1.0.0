import numpy as np
import sys
import hgf
import re
import os
import shutil

####################################################################
### PROBLEM SETUP, USER DEFINES INPUT/OUTPUT FOLDERS, VISCOSITY, ###
### NUMBER OF OMP THREADS, BOUNDARY CONDITIONS, GRID DIMENSIONS, ###
### AND SETS P-LEVEL FOR ILU PRECONDITIONER ########################
####################################################################

# FOLDER CONTAINING GRIDS
infolder = '/path/to/in/'
# DESTINATION FOR OUTPUTS
outfolder = '/path/to/out/'
# GRID DIMENSIONS
L = 1.
W = 1.
H = 1.

# PRINCIPAL FLOW DIRECTION, 0 - X, 1 - Y, 2 - Z
direction = 0

# SET VISCOSITY
visc = 0.001

# NUMBER OF OMP THREADS FOR USE IN PARALUTION LINEAR ALGEBRA
nThreads = 1

# SET ILU PRECONDITIONER LEVEL
prec = 6

####################################################################
### SETUP AND SWIG TRANSLATION, USER SHOULD NOT EDIT BELOW HERE ####
####################################################################

gridCount = 0
nGrids = len([name for name in os.listdir(infolder) if os.path.isfile(os.path.join(infolder, name))])

## Cycle through grids in directory ##
for root, dirs, filenames in os.walk(infolder):
  for f in filenames:
    gridCount = gridCount + 1
    gridF = open(os.path.join(root, f), 'r')
    gridF.seek(0)
    geoData = gridF.readline()
    gridin1 = np.genfromtxt(gridF, dtype = int)
    gridF.close()
    nxyz = re.findall(r'\d+', geoData)
    nx = int(nxyz[0])
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
                         L, W, H, direction, visc, nThreads, prec, nGrids, gridCount )

    ### Grab computed K ###
    KLoc = outfolder + 'Ks.dat'
    KFile = open('KConstantX.dat', 'r')
    KFile.seek(0)
    KStr = KFile.readline()
    KFlo = float(KStr)
    KFile.close()
    KOut = open(KLoc, 'a+')
    KOut.write('%f\t' % (KFlo))
    KOut.write(f)
    KOut.write('\n')
    KOut.close()

    ### Move flow solution file so not overwritten ###
    tarend = 'SOL' + f + '.dat'
    tarloc = outfolder + tarend
    shutil.copy('flowrun.dat', tarloc)
