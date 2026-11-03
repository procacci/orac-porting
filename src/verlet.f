      SUBROUTINE verlet(mass,ntap,xp0,yp0,zp0,vpx,vpy,vpz,fpxn,fpyn,fpzn
     &     ,dt)

************************************************************************
*                                                                      *
*  VERLET advance coordinates for one time step and velocities for     *
*  half time step. Coordinates and velocities arrays are overwritten.  *
*  This program is NOT a black box. Velocities are passed in common    *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      use omp_integr, only:omp_timing,time_in,tcpu_in
      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      integer      ntap
      real*8       dt,mass(ntap)
      REAL*8       xp0(ntap),yp0(ntap),zp0(ntap)
      REAL*8       vpx(ntap),vpy(ntap),vpz(ntap)
      REAL*8       fpxn(ntap),fpyn(ntap),fpzn(ntap)

!-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i
      REAL*8       ts2,tfact
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

!==================== EXECUTABLE STATEMENTS ============================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
      ts2=dt*dt

!$OMP DO SCHEDULE(STATIC) 
      do i=1,ntap
         xp0(i)=xp0(i)+dt*vpx(i)+fpxn(i)*0.5*ts2/mass(i)
         yp0(i)=yp0(i)+dt*vpy(i)+fpyn(i)*0.5*ts2/mass(i)
         zp0(i)=zp0(i)+dt*vpz(i)+fpzn(i)*0.5*ts2/mass(i)
         vpx(i)=vpx(i) + fpxn(i)*0.5*dt/mass(i)
         vpy(i)=vpy(i) + fpyn(i)*0.5*dt/mass(i)
         vpz(i)=vpz(i) + fpzn(i)*0.5*dt/mass(i)
      end do   
!OMP  ENDDO

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
      RETURN
      END

