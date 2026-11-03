      SUBROUTINE write_bonds(ktopol,fstep,top_bonds,lbnd,lbond,xp0,yp0
     &     ,zp0)

************************************************************************
*   Time-stamp: <97/07/08 19:43:05 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER top_bonds(*),lbnd(2,*),lbond,ktopol
      REAL*8  xp0(*),yp0(*),zp0(*),fstep

      REAL*8  bond(ntopol)
      COMMON /rag1/ bond
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i1,i,la,lb,n,j
      REAL*8  xr1,yr1,zr1,xr2,yr2,zr2,x21,y21,z21

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i1=1,top_bonds(1)
         i=top_bonds(1+i1)
         la=lbnd(1,i)
         lb=lbnd(2,i)
         xr1=xp0(la)
         yr1=yp0(la)
         zr1=zp0(la)
         xr2=xp0(lb)
         yr2=yp0(lb)
         zr2=zp0(lb)
         
         x21=xr2-xr1
         y21=yr2-yr1
         z21=zr2-zr1
         bond(i1)=DSQRT(x21**2+y21**2+z21**2)
      END DO
      DO i1=1,top_bonds(1),4
         n=3
         IF(top_bonds(1)-i1 .LT. 3) n=top_bonds(1)-i1
         WRITE(ktopol,1000) 'T',fstep,' Bonds ',(bond(j),top_bonds(1+j)
     &        ,j=i1,i1+n)
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1000  FORMAT(a1,f11.2,a7,1x,4(f12.5,2x,i5))
      RETURN
      END
