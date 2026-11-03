      SUBROUTINE set_const_co(co,dss,cnst)

************************************************************************
*   Time-stamp: <97/03/10 13:37:38 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 10 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER cnst(2,*)
      REAL*8  co(3,3),dss(*)

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      dss(1)=co(2,2)/co(1,1)
      dss(2)=co(3,3)/co(1,1)
      dss(3)=co(1,2)/co(1,1)
      dss(4)=co(1,3)/co(1,1)
      dss(5)=co(2,3)/co(1,1)
      cnst(1,1)=2
      cnst(2,1)=2
      cnst(1,2)=3
      cnst(2,2)=3
      cnst(1,3)=1
      cnst(2,3)=2
      cnst(1,4)=1
      cnst(2,4)=3
      cnst(1,5)=2
      cnst(2,5)=3

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
