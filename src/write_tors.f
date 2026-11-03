      SUBROUTINE write_tors(char,ktopol,fstep,top_ptors,ltors,ltorsion
     &     ,xp0,yp0,zp0)

************************************************************************
*   Time-stamp: <97/08/22 16:39:07 marchi>                             *
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

      INTEGER top_ptors(*),ltors(4,*),ltorsion,ktopol
      REAL*8  xp0(*),yp0(*),zp0(*),fstep
      CHARACTER*1 char

      REAL*8  tors(ntopol)
      COMMON /rag1/ tors
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i1,i,la,lb,lc,ld,n,j
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x12,x23
     &     ,x34,y12,y23,y34,z12,z23,z34,rsq12,rsq23,rsq34,rsp12,rsp23
     &     ,rsp34,pi,e12x,e12y,e12z,e23x,e23y,e23z,e34x,e34y,e34z,u1223x
     &     ,u1223y,u1223z,u2334x,u2334y,u2334z
      REAL*8  cb2,cb3,sb2,sb3,coa,aux,tsign

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      pi=4.0D0*DATAN(1.0D0)
      DO i1=1,top_ptors(1)
         i=top_ptors(1+i1)
         la=ltors(1,i)
         lb=ltors(2,i)
         lc=ltors(3,i)
         ld=ltors(4,i)
         xr1=xp0(la)
         yr1=yp0(la)
         zr1=zp0(la)
         xr2=xp0(lb)
         yr2=yp0(lb)
         zr2=zp0(lb)
         xr3=xp0(lc)
         yr3=yp0(lc)
         zr3=zp0(lc)
         xr4=xp0(ld)
         yr4=yp0(ld)
         zr4=zp0(ld)
         x12=xr1-xr2
         y12=yr1-yr2
         z12=zr1-zr2
         x23=xr2-xr3
         y23=yr2-yr3
         z23=zr2-zr3
         x34=xr3-xr4
         y34=yr3-yr4
         z34=zr3-zr4
         rsq12=x12**2+y12**2+z12**2
         rsq23=x23**2+y23**2+z23**2
         rsq34=x34**2+y34**2+z34**2
         rsp12=DSQRT(rsq12)
         rsp23=DSQRT(rsq23)
         rsp34=DSQRT(rsq34)
         e12x=x12/rsp12
         e12y=y12/rsp12
         e12z=z12/rsp12
         e23x=x23/rsp23
         e23y=y23/rsp23
         e23z=z23/rsp23
         e34x=x34/rsp34
         e34y=y34/rsp34
         e34z=z34/rsp34
         u1223x=e12y*e23z-e12z*e23y
         u1223y=e12z*e23x-e12x*e23z
         u1223z=e12x*e23y-e12y*e23x
         u2334x=e23y*e34z-e23z*e34y
         u2334y=e23z*e34x-e23x*e34z
         u2334z=e23x*e34y-e23y*e34x

         cb2=(x12*x23+y12*y23+z12*z23)/(rsp12*rsp23)
         cb3=(x34*x23+y34*y23+z34*z23)/(rsp34*rsp23)

         sb2=DSQRT(1.0d0-cb2**2)
         sb3=DSQRT(1.0d0-cb3**2)
         coa=(u1223x*u2334x+u1223y*u2334y+u1223z*u2334z)/sb2/sb3
         aux=(e12x*u2334x+e12y*u2334y+e12z*u2334z)/sb3
         tsign=DSIGN(1.0D0,aux)
         IF(coa.GT.1.0d0.AND.coa.LT.1.000001) coa=1.0d0
         IF(coa.LT.-1.0d0.AND.coa.GT.-1.000001) coa=-1.0d0
         coa=DACOS(coa)*180.0D0/pi
         tors(i1)=-(tsign-1.0D0)*180.0D0+tsign*coa
      END DO
      DO i1=1,top_ptors(1),4
         n=3
         IF(top_ptors(1)-i1 .LT. 3) n=top_ptors(1)-i1
         WRITE(ktopol,1000) 'T',fstep,' ',char,'-Tors ',(tors(j)
     &        ,top_ptors(1+j),j=i1,i1+n)
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1000  FORMAT(a1,f11.2,a1,a1,a6,4(f12.5,2x,i5))
      RETURN
      END
