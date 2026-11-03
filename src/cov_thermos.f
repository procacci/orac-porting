      SUBROUTINE cov_thermos(slv_exist,slt_exist,qmass,ntot,t)

************************************************************************
*   Time-stamp: <98/02/12 00:25:45 marchi>                             *
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

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,ntot(*)
      REAL*8  qmass(*),t
      LOGICAL slv_exist,slt_exist

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  c,omega

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      c=2.0D0*pi*2.997925D10

      omega=qmass(1)*c
      qmass(1)=2.0D0*DFLOAT(ntot(1))*boltz*t/omega**2
      qmass(1)=qmass(1)/(unitm*unitl**2)

      IF(slt_exist) THEN
         omega=qmass(2)*c
         qmass(2)=2.0D0*DFLOAT(ntot(2))*boltz*t/omega**2
         qmass(2)=qmass(2)/(unitm*unitl**2)
      ELSE
         qmass(2)=0.0D0
      END IF
      IF(slv_exist) THEN
         omega=qmass(3)*c
         qmass(3)=2.0D0*DFLOAT(ntot(3))*boltz*t/omega**2
         qmass(3)=qmass(3)/(unitm*unitl**2)
      ELSE
         qmass(3)=0.0D0
      END IF

      DO i=1,3
         IF(ntot(i) .EQ. 1) THEN
            qmass(i)=0.0D0
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
