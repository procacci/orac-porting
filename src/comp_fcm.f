      SUBROUTINE comp_fcm(nprot,protl,fpx,fpy,fpz,fcx,fcy,fcz,mass
     &     ,tmass,oc)

************************************************************************
*   Time-stamp: <97/06/26 13:04:10 marchi>                             *
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

      INTEGER nprot,protl(*)
      REAL*8  fpx(*),fpy(*),fpz(*),fcx(*),fcy(*),fcz(*),mass(*)
     &     ,tmass(*),oc(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,count,m,i1
      REAL*8  tt,xc,yc,zc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO i=1,nprot
         xc=0.0D0
         yc=0.0D0
         zc=0.0D0
         tt=1.0D0/tmass(i)
         m=protl(count+1)
         DO i1=1,m
            j=protl(count+1+i1)
            xc=xc+fpx(j)
            yc=yc+fpy(j)
            zc=zc+fpz(j)
         END DO
         fcx(i)=oc(1,1)*xc+oc(1,2)*yc+oc(1,3)*zc
         fcy(i)=oc(2,1)*xc+oc(2,2)*yc+oc(2,3)*zc
         fcz(i)=oc(3,1)*xc+oc(3,2)*yc+oc(3,3)*zc
         xc=xc*tt
         yc=yc*tt
         zc=zc*tt
         DO i1=1,m            
            j=protl(count+1+i1)
            fpx(j)=fpx(j)-mass(j)*xc
            fpy(j)=fpy(j)-mass(j)*yc
            fpz(j)=fpz(j)-mass(j)*zc
         END DO
         count=count+m+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
