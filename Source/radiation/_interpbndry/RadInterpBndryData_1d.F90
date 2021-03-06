#include "AMReX_BC_TYPES.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "RADINTERPBNDRYDATA_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 1
      
! ---------------------------------------------------------------
! ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPXLO(bdry,DIMS(bdry), &
                            lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                            mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(crse)
  integer  DIMDEC(cb)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(1)
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)

  integer  i, ic, n

  ic   = ARG_L1(cb)-1
  i    = lo(1)-1

  do n = 1, nvar
     ! interpolate to fine grid
     bdry(i,n) = crse(ic,n)
  enddo

  return
end subroutine FORT_BDINTERPXLO


! ---------------------------------------------------------------
! ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry), &
                             lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                             mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(cb)
  integer  DIMDEC(crse)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(1)
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)
  
  integer  i, ic, n

  ic   = ARG_H1(cb)+1
  i    = hi(1)+1
      
  do n = 1, nvar
     ! interpolate to fine grid
     bdry(i,n) = crse(ic,n)
  enddo

  return
end subroutine FORT_BDINTERPXHI
