      SUBROUTINE comp_cell_center(nato,xpc,ypc,zpc,xp0,yp0,zp0)

************************************************************************
*   Time-stamp: <97/12/05 19:01:35 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Dec  5 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato
      REAL*8 xp0(*),yp0(*),zp0(*),xpc,ypc,zpc

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nato
         xpc=xpc+xp0(i)
         ypc=ypc+yp0(i)
         zpc=zpc+zp0(i)
      END DO
      xpc=-xpc/DFLOAT(nato)
      ypc=-ypc/DFLOAT(nato)
      zpc=-zpc/DFLOAT(nato)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
