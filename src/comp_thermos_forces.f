      SUBROUTINE comp_thermos_forces(cpress,ndf_thermos,nato,ss_index
     &     ,nprot,co,vpx,vpy,vpz,vcmx,vcmy,vcmz,vco,mass,masscm,masspr,t
     &     ,fth)

************************************************************************
*   Time-stamp: <98/03/17 18:23:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Apr  5 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use omp_integr, only:omp_timing,time_in,tcpu_in
      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ndf_thermos(*),nato,nprot,ss_index(*)
      REAL*8  masspr,t,vpx(*),vpy(*),vpz(*),vcmx(*),vcmy(*),vcmz(*)
     &     ,mass(*),masscm(*),fth(*),vco(3,3),co(3,3)
      LOGICAL cpress

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,k,j
      REAL*8  tto(3),xd,yd,zd
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      if(omp_timing)  call timer(treal,tcpu1,tela1)
      tto(1)=0.0D0
      tto(2)=0.0D0
      tto(3)=0.0D0
!$OMP PARALLEL DO PRIVATE(xd,yd,zd) REDUCTION(+:tto(1))
      DO i=1,nprot
         xd=co(1,1)*vcmx(i)+co(1,2)*vcmy(i)+co(1,3)*vcmz(i)
         yd=                co(2,2)*vcmy(i)+co(2,3)*vcmz(i)
         zd=                                co(3,3)*vcmz(i)
         tto(1)=tto(1)+masscm(i)*(xd**2+yd**2+zd**2)
      END DO
!$OMP END PARALLEL DO

      IF(cpress) THEN
         DO i=1,3
            DO j=i,3
               tto(1)=tto(1)+masspr*vco(i,j)**2
            END DO
         END DO
      END IF

!$OMP PARALLEL DO PRIVATE(k) REDUCTION(+:tto)
      DO i=1,nato
         k=ss_index(i)+1
         tto(k)=tto(k)+mass(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
      END DO
!$OMP END PARALLEL DO
      
      DO i=1,3
        fth(i)=tto(i)-DFLOAT(ndf_thermos(i))*boltz*t/unite
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
      RETURN
      END
