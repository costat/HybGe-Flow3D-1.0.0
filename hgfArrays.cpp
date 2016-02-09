#include <iostream>
#include <cstdlib>
#include <vector>

#include "hgfMesh.hpp"
#include "hgfArrays.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

/* diffusionArrays performs a finite volume discretization
   for diffusion of a scalar variable. This is used
   to compute, component wise, discretizations of the diffusive
   term of the velocity in the Stokes
   or Navier-Stokes equations. */
void
diffusionArrays( const FluidMesh& Mesh, std::vector<int>& matIs, \
                 std::vector<int>& matJs, \
                 std::vector<double>& matVals, \
                 const std::vector<unsigned long>& velint, \
                 const std::vector<unsigned long>& fconn, \
                 const std::vector<double>& dxyz, \
                 double visc, int velShift )
{
  int cl, rcl;
  switch ( Mesh.DIM )
  {
    case 2 :
    {
      int nbrfaces[4];
      double dx[4], dy[4], val[5];
      /* Loop over interior nodes, ignores boundary, compute
         values and locations of nonzero entries along matrix row
         corresponding to each interior node. These are then stored
         in COO sparse matrix format */
      for (unsigned long arrayIndex = 0; arrayIndex < velint.size(); arrayIndex++) {
        cl = velint[arrayIndex];
        rcl = cl + velShift;
        for (int i = 0; i < 4; i++) {
          nbrfaces[i] = fconn[ idx2(cl, i, 4) ] - 1;
          dx[i] = dxyz[ idx2(nbrfaces[i], 0, 2) ];
          dy[i] = dxyz[ idx2(nbrfaces[i], 1, 2) ];
        }
        /* Values to be inserted in this row */
        val[0] = -visc * dxyz[ idx2(cl,0,2) ] / (0.5 * (dy[0] \
                                       + dxyz[ idx2(cl,1,2) ]));
        val[1] = -visc * dxyz[ idx2(cl,1,2) ] / (0.5 * (dx[1] \
                                       + dxyz[ idx2(cl,0,2) ]));
        val[2] = -visc * dxyz[ idx2(cl,0,2) ] / (0.5 * (dy[2] \
                                       + dxyz[ idx2(cl,1,2) ]));
        val[3] = -visc * dxyz[ idx2(cl,1,2) ] / (0.5 * (dx[3] \
                                       + dxyz[ idx2(cl,0,2) ]));
        val[4] = -(val[0] + val[1] + val[2] + val[3]);

        /* Insert values into assembly vector matVals */
        for (int i = 0; i < 4; i++)
        {
          matIs.push_back(rcl);
          matJs.push_back((nbrfaces[i] + velShift));
          matVals.push_back(val[i]);
        }
        matIs.push_back(rcl);
        matJs.push_back((cl + velShift));
        matVals.push_back(val[4]);
      }
      break;
    }
    case 3 :
    {
      int nbrfaces[6];
      double dx[6], dy[6], dz[6], val[7];
      /* Loop over interior nodes, ignores boundary, compute
         values and locations of nonzero entries along matrix row
         corresponding to each interior node. These are then stored
         in COO sparse matrix format */
      for (unsigned long arrayIndex = 0; arrayIndex < velint.size(); arrayIndex++) {
        cl = velint[arrayIndex];
        rcl = cl + velShift;
        for (int i = 0; i < 6; i++) {
          nbrfaces[i] = fconn[ idx2(cl, i, 6) ] - 1;
          dx[i] = dxyz[ idx2(nbrfaces[i], 0, 3) ];
          dy[i] = dxyz[ idx2(nbrfaces[i], 1, 3) ];
          dz[i] = dxyz[ idx2(nbrfaces[i], 2, 3) ];
        }
        /* Values to be inserted in this row */
        val[0] = -visc * dxyz[ idx2(cl,0,3) ] \
                       * dxyz[ idx2(cl,1,3) ] / (0.5 * (dz[0] \
                       + dxyz[ idx2(cl,2,3) ] ));
        val[1] = -visc * dxyz[ idx2(cl,1,3) ] \
                       * dxyz[ idx2(cl,2,3) ] / (0.5 * (dx[1] \
                       + dxyz[ idx2(cl,0,3) ] ));
        val[2] = -visc * dxyz[ idx2(cl,0,3) ] \
                       * dxyz[ idx2(cl,1,3) ] / (0.5 * (dz[2] \
                       + dxyz[ idx2(cl,2,3) ] ));
        val[3] = -visc * dxyz[ idx2(cl,1,3) ] \
                       * dxyz[ idx2(cl,2,3) ] / (0.5 * (dx[3] \
                       + dxyz[ idx2(cl,0,3) ] ));
        val[4] = -visc * dxyz[ idx2(cl,0,3) ] \
                       * dxyz[ idx2(cl,2,3) ] / (0.5 * (dy[4] \
                       + dxyz[ idx2(cl,1,3) ] ));
        val[5] = -visc * dxyz[ idx2(cl,0,3) ] \
                       * dxyz[ idx2(cl,2,3) ] / (0.5 * (dy[5] \
                       + dxyz[ idx2(cl,1,3) ] ));
        val[6] = -(val[0] + val[1] + val[2] + val[3] + val[4] + val[5]);

        /* Insert values into assembly vectors */
        for (int i = 0; i < 6; i++)
        {
          matIs.push_back(rcl);
          matJs.push_back((nbrfaces[i] + velShift));
          matVals.push_back(val[i]);
        }
        matIs.push_back(rcl);
        matJs.push_back((cl + velShift));
        matVals.push_back(val[6]);
      }
      break;
    }
  }
}

/* diffusionDrive organizes mesh information for each component
   call to diffusionArrays */
void
diffusionDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                 std::vector<int>& matJs, \
                 std::vector<double>& matVals, \
                 int Shift, double visc )
{
  int velShift;

  switch ( Shift )
  {
    case 0 :
    {
      const std::vector<unsigned long>& velint = Mesh.UInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.UFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.UCellWidths;
      velShift = 0;
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );
      break;
    }
    case 1 :
    {
      const std::vector<unsigned long>& velint = Mesh.VInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.VFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.VCellWidths;
      velShift = Mesh.DOF[1];
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );

      break;
    }
    case 2 :
    {
      const std::vector<unsigned long>& velint = Mesh.WInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.WFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.WCellWidths;
      velShift = Mesh.DOF[1] + Mesh.DOF[2];
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );

      break;
    }
  }
}

/* pressureGrad computes a finite volume discretization of the pressure
   gradient term in the Stokes and Navier-Stokes equations. */
void
pressureGrad ( const FluidMesh& Mesh, std::vector<int>& matIs, \
               std::vector<int>& matJs, std::vector<double>& matVals, \
               const std::vector<unsigned long>& velint, \
               const std::vector<unsigned long>& vcpn, \
               const std::vector<double>& dxyz, \
               int velShift, int done, int dtwo )
{

  int cl, rcl, totalShift, col[2];
  double val[2];

  switch ( Mesh.DIM )
  {
    case 3 :
    {
      totalShift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      /* Loop over interior nodes, ignores boundary, compute values
         and locations of nonzero entires along matrix row corresponding
         to each interior node. These are then stored in COO sparse matrix
         format */
      for (unsigned long arrayIndex = 0; arrayIndex < velint.size(); arrayIndex++) {
        cl = velint[arrayIndex];
        rcl = cl + velShift;
        val[0] = -dxyz[ idx2(cl, done, 3) ] * dxyz[ idx2(cl, dtwo, 3) ];
        val[1] = -val[0];
        col[0] = vcpn[ idx2( cl, 0, 2 ) ] - 1 + totalShift;
        col[1] = vcpn[ idx2( cl, 1, 2 ) ] - 1 + totalShift;
        // Insert values into assembly vector matVals
        for (int i = 0; i < 2; i++)
        {
          matIs.push_back(rcl);
          matJs.push_back(col[i]);
          matVals.push_back(val[i]);
        }
      }
      break;
    }
    case 2:
    {
      totalShift = Mesh.DOF[1] + Mesh.DOF[2];
      /* Loop over interior nodes, ignores boundary, compute values
         and locations of nonzero entires along matrix row corresponding
         to each interior node. These are then stored in COO sparse matrix
         format */
      for (unsigned long arrayIndex = 0; arrayIndex < velint.size(); arrayIndex++) {
        cl = velint[arrayIndex];
        rcl = cl + velShift;
        val[0] = -dxyz[ idx2(cl, done, 2) ];
        val[1] = -val[0];
        col[0] = vcpn[ idx2( cl, 0, 2 ) ] - 1 + totalShift;
        col[1] = vcpn[ idx2( cl, 1, 2 ) ] - 1 + totalShift;
        // Insert values into assembly vector matVals
        for (int i = 0; i < 2; i++)
        {
          matIs.push_back(rcl);
          matJs.push_back(col[i]);
          matVals.push_back(val[i]);
        }
      }
      break;
    }
  }
}

/* pressureGradDrive organizes mesh information for each component
   call to pressureGrad */
void
pressureGradDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                    std::vector<int>& matJs, \
                    std::vector<double>& matVals, \
                    int Shift )
{
  int velShift, done, dtwo;

  switch ( Shift )
  {
    case 0 :
    {
      const std::vector<unsigned long> &velint = Mesh.UInteriorCells;
      const std::vector<unsigned long> &vcpn = Mesh.UCellPressureNeighbor;
      const std::vector<double> &dxyz = Mesh.UCellWidths;
      velShift = 0;
      done = 1;
      dtwo = 2;
      pressureGrad ( Mesh, matIs, matJs, matVals, velint, vcpn, \
                     dxyz, velShift, done, dtwo );
      break;
    }
    case 1 :
    {
      const std::vector<unsigned long> &velint = Mesh.VInteriorCells;
      const std::vector<unsigned long> &vcpn = Mesh.VCellPressureNeighbor;
      const std::vector<double> &dxyz = Mesh.VCellWidths;
      velShift = Mesh.DOF[1];
      done = 0;
      dtwo = 2;
      pressureGrad ( Mesh, matIs, matJs, matVals, velint, vcpn, \
                     dxyz, velShift, done, dtwo );
      break;
    }
    case 2 :
    {
      const std::vector<unsigned long> &velint = Mesh.WInteriorCells;
      const std::vector<unsigned long> &vcpn = Mesh.WCellPressureNeighbor;
      const std::vector<double> &dxyz = Mesh.WCellWidths;
      velShift = Mesh.DOF[1] + Mesh.DOF[2];
      done = 0;
      dtwo = 1;
      pressureGrad ( Mesh, matIs, matJs, matVals, velint, vcpn, \
                             dxyz, velShift, done, dtwo );
      break;
    }
  }
}

/* continuityEq performs a finite volume discretization
   of the continuity equation */
void
continuityEq ( const FluidMesh& Mesh, std::vector<int>& matIs, \
               std::vector<int>& matJs, std::vector<double>& matVals, \
               const std::vector<unsigned long>& pcvn, \
               int velShift, int done, int dtwo )
{

  int totalShift, col[2];
  double val[2];

  switch ( Mesh.DIM )
  {
    case 3 :
    {
      totalShift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      /* Loop over all nodes, compute values
         and locations of nonzero entires along matrix row corresponding
         to each node. These are then stored in COO sparse matrix
         format */
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        val[0] = -1. * Mesh.PCellWidths[ idx2(cl, done, 3) ] \
                     * Mesh.PCellWidths[ idx2(cl, dtwo, 3) ];
        val[1] = -val[0];
        col[0] = pcvn[ idx2( cl, 0, 2 ) ] - 1 + velShift;
        col[1] = pcvn[ idx2( cl, 1, 2 ) ] - 1 + velShift;
        // Insert values into assembly vector
        for (int i = 0; i < 2; i++)
        {
          matVals.push_back(val[i]);
          matIs.push_back((totalShift + cl));
          matJs.push_back(col[i]);
        }
      }
      break;
    }
    case 2 :
    {
      totalShift = Mesh.DOF[1] + Mesh.DOF[2];
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        val[0] = -1. * Mesh.PCellWidths[ idx2(cl, done, 2) ];
        val[1] = -val[0];
        col[0] = pcvn[ idx2( cl, 0, 2 ) ] -1 + velShift;
        col[1] = pcvn[ idx2( cl, 1, 2 ) ] -1 + velShift;
        for (int i = 0; i < 2; i++)
        {
          matVals.push_back(val[i]);
          matIs.push_back((totalShift + cl));
          matJs.push_back(col[i]);
        }
      }
      break;
    }
  }
}

/* continuityDrive organizes mesh information for each component
   call to continuityEq */
void
continuityDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, \
                  std::vector<double>& matVals, \
                  int Shift )
{
  int velShift, done, dtwo;

  switch ( Shift )
  {
    case 0 :
    {
      const std::vector<unsigned long> &pcvn = Mesh.PressureCellUNeighbor;
      velShift = 0;
      done = 1;
      dtwo = 2;
      continuityEq ( Mesh, matIs, matJs, matVals, pcvn, velShift, \
                     done, dtwo );
      break;
    }
    case 1 :
    {
      const std::vector<unsigned long> &pcvn = Mesh.PressureCellVNeighbor;
      velShift = Mesh.DOF[1];
      done = 0;
      dtwo = 2;
      continuityEq ( Mesh, matIs, matJs, matVals, pcvn, velShift, \
                     done, dtwo );
      break;
    }
    case 2 :
    {
      const std::vector<unsigned long> &pcvn = Mesh.PressureCellWNeighbor;
      velShift = Mesh.DOF[1] + Mesh.DOF[2];
      done = 0;
      dtwo = 1;
      continuityEq ( Mesh, matIs, matJs, matVals, pcvn, velShift, \
                     done, dtwo );
      break;
    }
  }
}

/* StokesArray organizes the call to the previously defined functions
   to set up a finite volume discretization of the Stokes equations or
   the linear part of the Navier-Stokes equations in 2 or 3D.
   Boundary conditions for the velocity variables are not imposed, and
   the immersed boundary is not imposed. */
void
StokesArray ( const FluidMesh& Mesh, double visc, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals )
{
  switch ( Mesh.DIM ) {
    case 2 :
    {
      // Diffusion term in momentum equations
      diffusionDrive( Mesh, matIs, matJs, matVals, 0, visc );
      diffusionDrive( Mesh, matIs, matJs, matVals, 1, visc );

      // Pressure gradient in momentum equations
      pressureGradDrive( Mesh, matIs, matJs, matVals, 0 );
      pressureGradDrive( Mesh, matIs, matJs, matVals, 1 );

      // Arrays for continuity equation
      continuityDrive( Mesh, matIs, matJs, matVals, 0 );
      continuityDrive( Mesh, matIs, matJs, matVals, 1);

      break;
    }
    case 3 :
    {
      // Diffusion term in momentum equations
      diffusionDrive( Mesh, matIs, matJs, matVals, 0, visc );
      diffusionDrive( Mesh, matIs, matJs, matVals, 1, visc );
      diffusionDrive( Mesh, matIs, matJs, matVals, 2, visc );

      // Pressure gradient in momentum equations
      pressureGradDrive( Mesh, matIs, matJs, matVals, 0 );
      pressureGradDrive( Mesh, matIs, matJs, matVals, 1 );
      pressureGradDrive( Mesh, matIs, matJs, matVals, 2 );

      // Arrays for continuity equation
      continuityDrive( Mesh, matIs, matJs, matVals, 0 );
      continuityDrive( Mesh, matIs, matJs, matVals, 1 );
      continuityDrive( Mesh, matIs, matJs, matVals, 2 );

      break;
    }
  }
}
void
VelocityArray ( const FluidMesh& Mesh, double visc, std::vector<int>& matIs, \
                std::vector<int>& matJs, std::vector<double>& matVals, int component )
{
  int velShift;

  switch ( component )
  {
    case 0 :
    {
      const std::vector<unsigned long>& velint = Mesh.UInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.UFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.UCellWidths;
      velShift = 0;
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );
      break;
    }
    case 1 :
    {
      const std::vector<unsigned long>& velint = Mesh.VInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.VFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.VCellWidths;
      velShift = 0;
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );

      break;
    }
    case 2 :
    {
      const std::vector<unsigned long>& velint = Mesh.WInteriorCells;
      const std::vector<unsigned long>& fconn = Mesh.WFaceConnectivity;
      const std::vector<double>& dxyz = Mesh.WCellWidths;
      velShift = 0;
      diffusionArrays( Mesh, matIs, matJs, matVals, velint, fconn, \
                       dxyz, visc, velShift );

      break;
    }
  }
}
void
PoreNetworkArray( const PoreNetwork& pn, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  const std::vector<double>& Ks )
{
  double val[ 2*np.DIM + 1 ];
  int colId[ 2*np.DIM + 1 ];
  int pore;
  switch (np.DIM)
  {
    case 2 :
    {
      for (int ip = 0; ip < np.InteriorPores.size(); ip++)
      {
        pore = np.InteriorPores[ ip ];
        for (int dir = 0; dir < (np.DIM * 2); dir++) {
          colId[dir] = Throats[ idx2( pore, dir, (2*np.DIM) ) ];
        }
        colId[4] = pore;
        val[0] = -0.5 * ( Ks[ idx2( colId[0], 1, np.DIM ) ] + Ks[ idx2( pore, 1, np.DIM ) ] );
        val[1] = -0.5 * ( Ks[ idx2( colId[1], 0, np.DIM ) ] + Ks[ idx2( pore, 0, np.DIM ) ] );
        val[2] = -0.5 * ( Ks[ idx2( colId[2], 1, np.DIM ) ] + Ks[ idx2( pore, 1, np.DIM ) ] );
        val[3] = -0.5 * ( Ks[ idx2( colId[3], 0, np.DIM ) ] + Ks[ idx2( pore, 0, np.DIM ) ] );
        val[4] = -(val[0] + val[1] + val[2] + val[3]);
        for (int dir = 0; dir < (np.DIM * 2+1); dir++) {
          matIs.push_back( pore );
          matJs.push_back( colId[dir] );
          matVals.push_back( val[dir] );
        }
      }
      break;
    }
    default :
    {
      for (int ip = 0; ip < np.InteriorPores.size(); ip++)
      {
        pore = np.InteriorPores[ ip ];
        for (int dir = 0; dir < (np.DIM * 2); dir++) {
          colId[dir] = Throats[ idx2( pore, dir, (2*np.DIM) ) ];
        }
        colId[6] = pore;
        val[0] = -0.5 * ( Ks[ idx2( colId[0], 2, np.DIM ) ] + Ks[ idx2( pore, 2, np.DIM ) ] );
        val[1] = -0.5 * ( Ks[ idx2( colId[1], 0, np.DIM ) ] + Ks[ idx2( pore, 0, np.DIM ) ] );
        val[2] = -0.5 * ( Ks[ idx2( colId[2], 2, np.DIM ) ] + Ks[ idx2( pore, 2, np.DIM ) ] );
        val[3] = -0.5 * ( Ks[ idx2( colId[3], 0, np.DIM ) ] + Ks[ idx2( pore, 0, np.DIM ) ] );
        val[4] = -0.5 * ( Ks[ idx2( colId[4], 1, np.DIM ) ] + Ks[ idx2( pore, 1, np.DIM ) ] );
        val[5] = -0.5 * ( Ks[ idx2( colId[5], 1, np.DIM ) ] + Ks[ idx2( pore, 1, np.DIM ) ] );
        val[6] = -(val[0] + val[1] + val[2] + val[3] + val[4] + val[5]);
        for (int dir = 0; dir < (np.DIM * 2+1); dir++) {
          matIs.push_back( pore );
          matJs.push_back( colId[dir] );
          matVals.push_back( val[dir] );
        }
      }
      break;
    }
  }
}

