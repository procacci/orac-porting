      subroutine correc(vpx,vpy,vpz,fpx,fpy,fpz,mass,ntap,dt)

************************************************************************
*                                                                      *
*     CORREC advances or velocities for half time step or corrects     *
*     velocities al full time step.                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use omp_integr, only:omp_timing,time_in,tcpu_in
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8       fpx(ntap),fpy(ntap),fpz(ntap)
      REAL*8       vpx(ntap),vpy(ntap),vpz(ntap)
      integer      ntap
      real*8       dt,mass(ntap)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i,itask
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*==================== EXECUTABLE STATEMENTS ============================
      
*=======================================================================
*  advance atomic velocities
*=======================================================================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
!$OMP DO SCHEDULE(STATIC) 
      DO i=1,ntap
        vpx(i)=vpx(i) + fpx(i)*0.5*dt/mass(i)
        vpy(i)=vpy(i) + fpy(i)*0.5*dt/mass(i)
        vpz(i)=vpz(i) + fpz(i)*0.5*dt/mass(i)
      END DO
!OMP  ENDDO 
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
      RETURN
      END
