      SUBROUTINE tr_inbox(xp0,yp0,zp0,xp1,yp1,zp1,mass,nprot,protl)

************************************************************************
*   Time-stamp: <98/03/19 16:40:12 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul 20 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use cpropar, only:ascii_nocell
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  xp0(*),yp0(*),zp0(*),xp1(*),yp1(*),zp1(*),mass(*)
      INTEGER nprot,protl(*)

*------------------------- LOCAL VARIABLES ----------------------------*
      
      INTEGER i,j,m,n,count
      REAL*8  xcm,ycm,zcm,xmap,ymap,zmap,tmass,xcm1,ycm1,zcm1
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*     origin is in the com of 1st molecule.
      if(ascii_nocell) THEN
        count=0
        DO i=1,nprot
          m=protl(count+1)
          DO j=1,m
            n=protl(count+1+j)
            xp1(n)=xp0(n)
            yp1(n)=yp0(n)
            zp1(n)=zp0(n)
          END DO
          count=count+m+1
        END DO
        RETURN
      END IF
      m=protl(1)
      xcm=0.0D0
      ycm=0.0D0
      zcm=0.0D0
      tmass=0.0D0
      DO j=1,m
        n=protl(1+j)
        xcm=xcm+xp0(n)*mass(n)
        ycm=ycm+yp0(n)*mass(n)
        zcm=zcm+zp0(n)*mass(n)
        tmass=tmass+mass(n)
      END DO
      xcm1=xcm/tmass
      ycm1=ycm/tmass
      zcm1=zcm/tmass
*     now shift all cordinates to new origin 
      count=0
      DO i=1,nprot
         m=protl(count+1)
         DO j=1,m
           n=protl(count+1+j)
           xp1(n)=xp0(n)-xcm1
           yp1(n)=yp0(n)-ycm1
           zp1(n)=zp0(n)-zcm1
         END DO
         count=count+m+1
      END DO
      count=0
*     move all molecules in the BOX (first solute molecule at the center)
      DO i=1,nprot
         m=protl(count+1)
         xcm=0.0D0
         ycm=0.0D0
         zcm=0.0D0
         tmass=0.0D0
         DO j=1,m
            n=protl(count+1+j)
            xcm=xcm+xp1(n)*mass(n)
            ycm=ycm+yp1(n)*mass(n)
            zcm=zcm+zp1(n)*mass(n)
            tmass=tmass+mass(n)
         END DO
         xcm=xcm/tmass
         ycm=ycm/tmass
         zcm=zcm/tmass
         xmap=2.0D0*PBC(xcm)
         ymap=2.0D0*PBC(ycm)
         zmap=2.0D0*PBC(zcm)
         DO j=1,m
            n=protl(count+1+j)
            xp1(n)=xp1(n)-xmap
            yp1(n)=yp1(n)-ymap
            zp1(n)=zp1(n)-zmap
         END DO
         count=count+m+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
