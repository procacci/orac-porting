      SUBROUTINE coordinate_spline_init(time,fstep,length)

************************************************************************
*   Time-stamp: <97/11/28 14:19:04 marchi>                             *
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

      INTEGER length
      REAL*8  time(*),fstep

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,length
         time(i)=fstep*DFLOAT(i-1)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
