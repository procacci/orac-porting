      SUBROUTINE low_up(string,ndim)

************************************************************************
*   Time-stamp: <95/04/03 13:32:35 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr  3 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ndim
      CHARACTER*1 string(ndim)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,nstr

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=1,ndim
         nstr=ICHAR(string(n))
         IF(nstr .GE. 97 .AND. nstr .LE. 122) THEN
            string(n)=CHAR(nstr-32)
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
