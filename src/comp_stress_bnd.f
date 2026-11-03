      SUBROUTINE comp_stress_bnd(ntap,atomp,virial,co,xcm,ycm,zcm,fx,fy
     &     ,fz)

************************************************************************
*   Time-stamp: <98/01/16 21:36:25 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER ntap,atomp(*)
      REAL*8  xcm(*),ycm(*),zcm(*),fx(*),fy(*),fz(*)
      REAL*8  co(3,3),virial(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,p1
      REAL*8  xc,yc,zc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,ntap
         p1=atomp(i)
         xc=co(1,1)*xcm(p1)+co(1,2)*ycm(p1)+co(1,3)*zcm(p1)
         yc=co(2,1)*xcm(p1)+co(2,2)*ycm(p1)+co(2,3)*zcm(p1)
         zc=co(3,1)*xcm(p1)+co(3,2)*ycm(p1)+co(3,3)*zcm(p1)
         virial(1,1)=virial(1,1)+xc*fx(i)
         virial(1,2)=virial(1,2)+yc*fx(i)
         virial(1,3)=virial(1,3)+zc*fx(i)
         virial(2,1)=virial(2,1)+xc*fy(i)
         virial(2,2)=virial(2,2)+yc*fy(i)
         virial(2,3)=virial(2,3)+zc*fy(i)
         virial(3,1)=virial(3,1)+xc*fz(i)
         virial(3,2)=virial(3,2)+yc*fz(i)
         virial(3,3)=virial(3,3)+zc*fz(i)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
