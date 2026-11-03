      subroutine correc_gmgp(vpx,vpy,vpz,gmgp,nmol)

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

      REAL*8       gmgp(3,3)
      REAL*8       vpx(*),vpy(*),vpz(*)
      integer      nmol


*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i
      real*8       xc,yc,zc
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*  advance atomic velocities
*=======================================================================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
!OMP PARALLEL DO
      do i=1,nmol
         xc=gmgp(1,1)*vpx(i)+gmgp(1,2)*vpy(i)+gmgp(1,3)*vpz(i)
         yc=gmgp(2,1)*vpx(i)+gmgp(2,2)*vpy(i)+gmgp(2,3)*vpz(i)
         zc=gmgp(3,1)*vpx(i)+gmgp(3,2)*vpy(i)+gmgp(3,3)*vpz(i)
         vpx(i)=xc
         vpy(i)=yc
         vpz(i)=zc
      end do  
!OMP END PARALLEL DO 
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF

      RETURN
      END
