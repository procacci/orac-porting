      SUBROUTINE change_origin(dir,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xcm
     &     ,ycm,zcm,co)

************************************************************************
*   Time-stamp: <97/06/26 13:01:55 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar  7 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER dir,nprot,protl(*)
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*),lx(*),ly(*),lz(*),xcm(*),ycm(
     &     *),zcm(*)

*------------------------- LOCAL VARIABLES ----------------------------*


      INTEGER i,j,count,m,i1
      REAL*8  xc,yc,zc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(dir .EQ. 1) THEN
         count=0
         DO i=1,nprot
            xc=co(1,1)*xcm(i)+co(1,2)*ycm(i)+co(1,3)*zcm(i)
            yc=co(2,1)*xcm(i)+co(2,2)*ycm(i)+co(2,3)*zcm(i)
            zc=co(3,1)*xcm(i)+co(3,2)*ycm(i)+co(3,3)*zcm(i)
            m=protl(1+count)
            DO i1=1,m
               j=protl(count+1+i1)
               lx(j)=xp0(j)-xc
               ly(j)=yp0(j)-yc
               lz(j)=zp0(j)-zc
            END DO
            count=count+m+1
         END DO
      ELSE IF(dir .EQ. -1) THEN
         count=0
         DO i=1,nprot
            xc=co(1,1)*xcm(i)+co(1,2)*ycm(i)+co(1,3)*zcm(i)
            yc=co(2,1)*xcm(i)+co(2,2)*ycm(i)+co(2,3)*zcm(i)
            zc=co(3,1)*xcm(i)+co(3,2)*ycm(i)+co(3,3)*zcm(i)
            m=protl(count+1)
            DO i1=1,m
               j=protl(count+1+i1)
               xp0(j)=lx(j)+xc
               yp0(j)=ly(j)+yc
               zp0(j)=lz(j)+zc
            END DO
            count=count+m+1
         END DO
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
