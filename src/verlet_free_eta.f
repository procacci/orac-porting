      SUBROUTINE verlet_free_eta(n,eta,etap,dt)

************************************************************************
*                                                                      *
*  VERLET advance coordinates for one time step and velocities for     *
*  half time step. Coordinates and velocities arrays are overwritten.  *
*  This program is NOT a black box. Velocities are passed in common    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      integer      n
      real*8       dt
      REAL*8       eta(*),etap(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i

*==================== EXECUTABLE STATEMENTS ============================

      do i=1,n
         eta(i)=eta(i)+dt*etap(i)
      end do   
      RETURN
      END

