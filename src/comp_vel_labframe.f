      SUBROUTINE comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co,nprot
     &     ,protl,massa,massb,vpcmx,vpcmy,vpcmz)

************************************************************************
*   Time-stamp: <98/01/21 15:14:00 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Dec  8 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  vpx(*),vpy(*),vpz(*),vpx1(*),vpy1(*),vpz1(*),vpcmx(*)
     &     ,vpcmy(*),vpcmz(*),massa(*),massb(*),co(3,3)
      INTEGER nprot,protl(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,count,m,i1
      REAL*8  mtot,xc,yc,zc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO j=1,nprot
         m=protl(1+count)
         xc=co(1,1)*vpcmx(j)+co(1,2)*vpcmy(j)+co(1,3)*vpcmz(j)
         yc=co(2,1)*vpcmx(j)+co(2,2)*vpcmy(j)+co(2,3)*vpcmz(j)
         zc=co(3,1)*vpcmx(j)+co(3,2)*vpcmy(j)+co(3,3)*vpcmz(j)
         DO i=1,m
            i1=protl(count+1+i)
            vpx1(i1)=vpx(i1)+xc
            vpy1(i1)=vpy(i1)+yc
            vpz1(i1)=vpz(i1)+zc
         END DO
         count=count+m+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
