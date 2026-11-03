      SUBROUTINE zero0(a,n)

************************************************************************
*   Time-stamp: <97/11/28 14:19:20 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Nov 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      REAL*8  a(*)
      INTEGER n

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,n
         a(i)=0.0D0
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
