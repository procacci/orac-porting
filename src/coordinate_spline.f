      SUBROUTINE coordinate_spline(spline_x,spline_y,xpc,lengtha,divide)

************************************************************************
*   Time-stamp: <97/12/05 17:22:35 marchi>                             *
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

      REAL*8  spline_x(*),spline_y(4,*),xpc(0:*)
      INTEGER lengtha,divide

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,i1,ibcbeg,ibcend,length

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      length=lengtha/divide
      DO i=0,length-1
         i1=i+1
         spline_y(1,i1)=xpc(i)
      END DO
      ibcbeg=0
      ibcend=0
      CALL cubspl(spline_x,spline_y,length,ibcbeg,ibcend)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
