      SUBROUTINE correc_scale_co(n,v,etap,dt)

************************************************************************
*   Time-stamp: <97/04/07 12:25:05 marchi>                             *
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
      REAL*8  etap,dt,v(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  scale

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      scale=DEXP(-dt*etap)
      DO i=1,n
         v(i)=v(i)*scale
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
