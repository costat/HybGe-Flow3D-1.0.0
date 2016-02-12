import numpy as np
import sys
import hgf
import re

##############################################################################
### PROBLEM SETUP, USER DEFINES GRID, VISCOSITY, AND NUMBER OF OMP THREADS ###
### AND P-LEVEL FOR PRECONDITIONER ###########################################
##############################################################################

# GRID INFORMATION. USER PROVIDES PATH .DAT FILE CONTAINING
# VOXEL ARRAY OF 0S 1S AND 2S.
# ALSO, USER PROVIDES TOTAL GRID LENGTHS IN EACH DIRECTION.
gridfiles = './grids/2dsquare.dat'
L = 1.
W = 1.
H = 1.

# PRINCIPAL FLOW DIRECTION, 0 - X, 1 - Y, 2 - Z, SINGLE FLOW DIRECTION SOLVES,
# PRODUCES CONSTANT K FOR USE IN PORE-NETWORK THROATS
# 3 - ALL DIRECTIONS, PRODUCES UPSCALED K TENSOR
direction = 4

# SET VISCOSITY
visc = 1

# NUMBER OF OMP THREADS FOR USE IN PARALUTION LINEAR ALGEBRA
nThreads = 4

# SET SOLVER PARAMETERS: ILU PRECONDITIONER LEVEL,
# ABSOLUTE AND RELATIVE RESIDUAL TOLERANCES, AND MAXIMUM ITERATIONS
prec = 0
tolAbs = 1e-8;
tolRel = 1e-8;
maxIt = 3000;
solver = 0;
relax = 0;
MX = 2;
MY = 2
MZ = 0;

##########################################################
### SWIG TRANSLATION, USER SHOULD NOT EDIT BELOW HERE ####
##########################################################

gridF = open( gridfiles, 'r' )
gridF.seek(0) # make sure we read the file from the beginning
geoData = gridF.readline() # read the first line containing nx, ny, nz info
gridin1 = np.genfromtxt(gridF, dtype = int) # import the rest as a numpy array
gridF.close()
nxyz = re.findall(r'\d+', geoData)
nx = int(nxyz[0])
ny = int(nxyz[1])
nz = int(nxyz[2])

if nz:
    gridin = np.zeros((nx, ny, nz), dtype = int)
    for zLevel in range(0,nz):
        gridin[:,:,zLevel] = gridin1[zLevel*nx:(zLevel+1)*nx][:]
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

hgf.hgfDrive ( gridin, gridin_ldi2, gridin_ldi3, nx, ny, nz, \
               L, W, H, direction, visc, nThreads, prec, 1, 1, \
               tolAbs, tolRel, maxIt, MX, MY, MZ, solver, relax )
