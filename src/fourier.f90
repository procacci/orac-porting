#include "config.H"
module fourier
  implicit none
  
  !============================================================================
  !---- Parameters of the sum over the reciprocal lattice ---------------------
  !============================================================================
  
  INTEGER :: lmax,mmax,nmax
  INTEGER :: lofset,mofset,nofset,lmaxt,mmaxt,nmaxt,mmint,nmint
  
  PARAMETER ( &
       LMAX = _LMAX_, MMAX = _MMAX_, NMAX = _NMAX_, &
       lofset = 1, mofset = MMAX + 1, nofset = NMAX + 1, &
       lmaxt = LMAX + 1, mmaxt = 2*MMAX + 1,  &
       nmaxt = 2*NMAX + 1)
  
  !============================================================================
  
end module fourier
