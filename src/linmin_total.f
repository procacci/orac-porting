      SUBROUTINE linmin_total(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,xp0,yp0,zp0,xp1,yp1,zp1,fpx,fpy,fpz,fpx1,fpy1,fpz1,n
     &     ,fret)

************************************************************************
*   Time-stamp: <98/03/06 18:20:36 marchi>                             *
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

      INTEGER mapnl(*),mapdn(2,*),nmapdn(*)
      INTEGER nnlpp0(*),nnlpp(*)
      REAL*8  fudgec
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xp1(*),yp1(*)
     &     ,zp1(*),fpx1(*),fpy1(*),fpz1(*),fret
      INTEGER n

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  brent_total,f1dim_total,ax,bx,xx,fa,fb,fc,xmin,tol,fx,x
      EXTERNAL brent_total,f1dim_total
      DATA tol/1.0D-4/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL dcopy(n,xp0,1,xp1,1)
      CALL dcopy(n,yp0,1,yp1,1)
      CALL dcopy(n,zp0,1,zp1,1)
      CALL dcopy(n,fpx,1,fpx1,1)
      CALL dcopy(n,fpy,1,fpy1,1)
      CALL dcopy(n,fpz,1,fpz1,1)

      ax=0.0D0
      xx=1.0D0
      CALL mnbrakf_total(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,ax,xx,bx,fa,fb,fc,f1dim_total,xp1,yp1,zp1,fpx1
     &     ,fpy1,fpz1)
      fret=brent_total(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,ax,xx,bx,f1dim_total,tol,xmin,xp1,yp1,zp1
     &     ,fpx1,fpy1,fpz1)

      CALL daxpy(n,xmin,fpx,1,xp0,1)
      CALL daxpy(n,xmin,fpy,1,yp0,1)
      CALL daxpy(n,xmin,fpz,1,zp0,1)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(' Energy = ',f13.5)
      RETURN
      END
