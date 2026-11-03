      SUBROUTINE correc_etap(n,etap,fth,qmass,dt)

************************************************************************
*   Time-stamp: <97/04/07 12:26:36 marchi>                             *
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

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n
      REAL*8  etap(*),qmass(*),fth(*),dt

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,n
         IF(DABS(qmass(i)) .GT. 1.0D-7) THEN
            etap(i)=etap(i)+dt*fth(i)/qmass(i)
         ELSE
            etap(i)=0.0D0
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
