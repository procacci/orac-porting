      SUBROUTINE correc_scale(noff,n,vx,vy,vz,etap,dt)

************************************************************************
*   Time-stamp: <97/04/07 12:26:09 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Apr  6 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use omp_integr, only:omp_timing,time_in,tcpu_in
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,noff
      REAL*8  etap,dt,vx(*),vy(*),vz(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  scale
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      if(omp_timing)  call timer(treal,tcpu1,tela1)
      scale=DEXP(-dt*etap)
!$OMP PARALLEL DO  
      DO i=noff,n
         vx(i)=vx(i)*scale
         vy(i)=vy(i)*scale
         vz(i)=vz(i)*scale
      END DO
!$OMP END PARALLEL DO

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
