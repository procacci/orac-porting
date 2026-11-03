      RECURSIVE SUBROUTINE perc_bondh(index,concth,m1,mask)

************************************************************************
*   Time-stamp: <97/09/12 12:44:46 marchi>                             *
*                                                                      *
*   Percolate through a bond network including only connections        *
*   with hydrogens                                                     *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER index,m1
      INTEGER concth(m1,*)
      LOGICAL mask(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,coordi,ia

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      mask(index)=.FALSE.
      coordi=concth(index,1)
      DO i=2,coordi+1
         ia=concth(index,i)
         IF(mask(ia)) THEN
            CALL perc_bondh(ia,concth,m1,mask)
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
