      SUBROUTINE comp_dip2(x0,y0,z0,charge,dip,nts)

************************************************************************
*   Time-stamp: <97/02/28 12:50:57 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Sep  6 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      REAL*8 dip(3),x0(*),y0(*),z0(*),charge(*)
      INTEGER nts

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 xc,yc,zc
      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      dip(1)=0.0D0
      dip(2)=0.0D0
      dip(3)=0.0D0
      DO i=1,nts
         xc=x0(i)
         yc=y0(i)
         zc=z0(i)
         dip(1)=dip(1)+charge(i)*xc
         dip(2)=dip(2)+charge(i)*yc
         dip(3)=dip(3)+charge(i)*zc
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
