      SUBROUTINE mts_plotp(beta,mback,nbone,xp0,yp0,zp0,ntap,res,ma,
     x                type)

************************************************************************
*                                                                      *
*     Print a snapshot of the backbone of the solute molecule.         *
*                                                                      *
*---- Last update 01/04/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley, 1991                      *
*                                                                      *
*     EXTERNALS NONE.                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,ma,nbone,mback(*)
      INTEGER res(ma,2)
      REAL*8  xp0(*),yp0(*),zp0(*)
      CHARACTER*7 beta(*)
      CHARACTER*8 type(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l,m,m1,m2,mb,mc,indx
      REAL*8  xb,yb,zb,charge
      CHARACTER*1 chia,hydr,carb,oxyg,nitr,sulph,phosp
      LOGICAL ok

*==================== EXECUTABLE STATEMENTS ============================

      charge=0.0D0
      k=0
      j=1
      write(kplot,*) nbone
      write (kplot,*) 'backbone'
      DO m=1,nbone
          i=mback(m)
          mb=res(i,2)
          mc=res(i,1)
          xb=xp0(i)
          yb=yp0(i)
          zb=zp0(i)
          WRITE(kplot,1) beta(i)(1:1),xb,yb,zb,charge
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

1     FORMAT(5x,a1,4f10.4)
      RETURN
      END

