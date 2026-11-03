      SUBROUTINE zero_gofr(maxint,krdf,ngrdon,offset)

************************************************************************
*   Time-stamp: <97/08/06 16:27:26 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Dec  7 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER maxint,ngrdon,offset,krdf(maxint,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ngrdon=0
      DO i=1,3+offset
         DO j=1,maxint
            krdf(j,i)=0
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
