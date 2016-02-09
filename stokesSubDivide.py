import numpy as np
import sys
import hgf
import re
import os
import shutil
import math

# PARENT GRID LOCATION
gridfiles = './grids/2dsquare.dat'
# GRID DIMENSIONS
L = 1.
W = 1.
H = 1.

# DOMAIN SUBDIVISION INFO
MX = 2
MY = 2
MZ = 0
totalGrids = MX * MY

# SET VISCOSITY
visc = 1.

# NUMBER OF OMP THREADS FOR USE IN PARALUTION LINEAR ALGEBRA
nThreads = 1

# SET SOLVER PARAMETERS: ILU PRECONDITIONER LEVEL,
# ABSOLUTE AND RELATIVE RESIDUAL TOLERANCES, AND MAXIMUM ITERATIONS
prec = 1
tolAbs = 1e-12;
tolRel = 1e-12;
maxIt = 5000;
solver = 0;
relax = 0;

###################################
### SETUP AND SWIG TRANSLATION ####
###################################

gridF = open( gridfiles, 'r' )
gridF.seek(0)
geoData = gridF.readline()
gridin1 = np.genfromtxt(gridF, dtype = int)
gridF.close()
nxyz = re.findall(r'\d+', geoData)
nx = int(nxyz[0])
ny = int(nxyz[1])
nz = int(nxyz[2])
subgrid = 0
gridCount = 0

if nz:
  print 'Not yet implemented for 3D grids.'

elif not nz:
  gridin1 = np.transpose(gridin1)
  yStart = 0
  nyRemainder = ny
  yCount = 0
  gridTotal = MX * MY * 2
  for yCount in range(0, MY):
    xStart = 0
    nxRemainder = nx
    xCount = 0
    nyG = int(round(float(nyRemainder)/(MY-yCount)))
    yLim = [yStart, yStart + nyG]
    for xCount in range(0, MX):
      gridCount = gridCount + 1
      subgrid = subgrid+1
      nxG = int(round(float(nxRemainder)/(MX-xCount)))
      xLim = [xStart, xStart + nxG]
      gridin = gridin1[xLim[0]:xLim[1],yLim[0]:yLim[1]]
      gridin_ldi1, gridin_ldi2 = gridin.shape
      gridin = np.array(gridin, dtype = np.uint64).ravel()
      gridin_ldi3 = 0
      lG = L * (float(nxG)/nx)
      wG = W * (float(nyG)/ny)

      hgf.hgfDrive( gridin, gridin_ldi2, gridin_ldi3, nxG, nyG, nz, \
                    lG, wG, H, 0, visc, nThreads, prec, gridTotal, gridCount, \
                    tolAbs, tolRel, maxIt, solver, relax )
      sol1 = 'SOL_grid%d_x.dat' % (subgrid)
      shutil.copy('flowrun.dat', sol1)

      Ksout = 'subGridKs.dat'
      KFile = open('KConstantX.dat', 'r')
      KFile.seek(0)
      KStr = KFile.readline()
      KFlo = float(KStr)
      KFile.close()
      KOut = open(Ksout, 'a+')
      KOut.write('%f\t' % (KFlo))
      KOut.write('Subgrid %d, x-flow' % (subgrid))
      KOut.write('\n')
      KOut.close()

      gridCount = gridCount + 1
      hgf.hgfDrive( gridin, gridin_ldi2, gridin_ldi3, nxG, nyG, nz, \
                    lG, wG, H, 1, visc, nThreads, prec, gridTotal, gridCount, \
                    tolAbs, tolRel, maxIt, solver, relax )
      sol2 = 'SOL_grid%d_y.dat' % (subgrid)
      shutil.copy('flowrun.dat', sol2)

      Ksout = 'subGridKs.dat'
      KFile = open('KConstantY.dat', 'r')
      KFile.seek(0)
      KStr = KFile.readline()
      KFlo = float(KStr)
      KFile.close()
      KOut = open(Ksout, 'a+')
      KOut.write('%f\t' % (KFlo))
      KOut.write('Subgrid %d, y-flow' % (subgrid))
      KOut.write('\n')
      KOut.close()

      xStart = xStart + nxG
      nxRemainder = nx - xStart
      xCount = xCount + 1
    yStart = yStart + nyG
    nyRemainder = ny - yStart
    yCount = yCount + 1
