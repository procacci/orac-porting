      SUBROUTINE covchg(charge,nato)

************************************************************************
*   Time-stamp: <95/03/27 19:34:41 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 27 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nato
      REAL*8  charge(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=1,nato
         charge(n)=charge(n)/DSQRT(unitc)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
