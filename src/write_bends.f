      SUBROUTINE write_bends(ktopol,fstep,top_bends,lbndg,lbend,xp0,yp0
     &     ,zp0)

************************************************************************
*   Time-stamp: <97/07/08 19:43:27 marchi>                             *
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

      INTEGER top_bends(*),lbndg(3,*),lbend,ktopol
      REAL*8  xp0(*),yp0(*),zp0(*),fstep

      REAL*8  bend(ntopol)
      COMMON /rag1/ bend
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i1,i,la,lb,lc,n,j
      REAL*8  xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,x12,y12,z12,x32,y32
     &     ,z32,rs12,rs32,dcc2,cb,pi

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      pi=4.0D0*DATAN(1.0D0)
      DO i1=1,top_bends(1)
         i=top_bends(1+i1)
         la=lbndg(1,i)
         lb=lbndg(2,i)
         lc=lbndg(3,i)
         xr1=xp0(la)
         yr1=yp0(la)
         zr1=zp0(la)
         xr2=xp0(lb)
         yr2=yp0(lb)
         zr2=zp0(lb)
         xr3=xp0(lc)
         yr3=yp0(lc)
         zr3=zp0(lc)

         x12=xr1-xr2
         y12=yr1-yr2
         z12=zr1-zr2
         x32=xr3-xr2
         y32=yr3-yr2
         z32=zr3-zr2

         rs12=x12**2+y12**2+z12**2
         rs32=x32**2+y32**2+z32**2
         dcc2=DSQRT(rs12*rs32)
         cb=(x12*x32+y12*y32+z12*z32)/dcc2
         bend(i1)=DACOS(cb)*180.0D0/pi
      END DO
      DO i1=1,top_bends(1),4
         n=3
         IF(top_bends(1)-i1 .LT. 3) n=top_bends(1)-i1
         WRITE(ktopol,1000) 'T',fstep,' Bends ',(bend(j),top_bends(1+j)
     &        ,j=i1,i1+n)
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1000  FORMAT(a1,f11.2,a7,1x,4(f12.5,2x,i5))
      RETURN
      END
