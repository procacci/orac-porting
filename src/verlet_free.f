      SUBROUTINE verlet_free(ntap,xp0,yp0,zp0,vpx,vpy,vpz,dt)

************************************************************************
*                                                                      *
*  VERLET advance coordinates for one time step and velocities for     *
*  half time step. Coordinates and velocities arrays are overwritten.  *
*  This program is NOT a black box. Velocities are passed in common    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use omp_integr, only:omp_timing,time_in,tcpu_in
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      integer      ntap
      real*8       dt
      REAL*8       xp0(*),yp0(*),zp0(*)
      REAL*8       vpx(*),vpy(*),vpz(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8       tfact
      INTEGER      i
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*==================== EXECUTABLE STATEMENTS ============================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
!$OMP PARALLEL DO
      do i=1,ntap
         xp0(i)=xp0(i)+dt*vpx(i)
         yp0(i)=yp0(i)+dt*vpy(i)
         zp0(i)=zp0(i)+dt*vpz(i)
      end do   
!$OMP END PARALLEL DO
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
      RETURN
      END

