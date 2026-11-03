      SUBROUTINE linmin_adjust_bonds(xp0,yp0,zp0,xp1,yp1,zp1,fpx,fpy,fpz
     &     ,fpx1,fpy1,fpz1,n,fret)

************************************************************************
*   Time-stamp: <97/02/23 12:55:07 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 23 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xp1(*),yp1(*)
     &     ,zp1(*),fpx1(*),fpy1(*),fpz1(*),fret
      INTEGER n

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  brent,f1dim_adjust_bonds,ax,bx,xx,fa,fb,fc,xmin,tol
      EXTERNAL brent,f1dim_adjust_bonds
      DATA tol/1.0D-4/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,n
         xp1(i)=xp0(i)
         yp1(i)=yp0(i)
         zp1(i)=zp0(i)
         fpx1(i)=fpx(i)
         fpy1(i)=fpy(i)
         fpz1(i)=fpz(i)
      END DO

      ax=0.0D0
      xx=1.0D0
      CALL mnbrakf(ax,xx,bx,fa,fb,fc,f1dim_adjust_bonds,xp1,yp1,zp1
     &     ,fpx1,fpy1,fpz1)

      fret=brent(ax,xx,bx,f1dim_adjust_bonds,tol,xmin,xp1,yp1,zp1,fpx1
     &     ,fpy1,fpz1)

      DO i=1,n
         fpx(i)=xmin*fpx(i)
         fpy(i)=xmin*fpy(i)
         fpz(i)=xmin*fpz(i)
         xp0(i)=xp0(i)+fpx(i)
         yp0(i)=yp0(i)+fpy(i)
         zp0(i)=zp0(i)+fpz(i)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
