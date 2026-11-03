      SUBROUTINE write_gyr(kgyr,fstep,protl,nprot,ss_index,xp0,yp0,zp0)

************************************************************************
*   Time-stamp: <97/07/15 23:12:52 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul 15 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kgyr,nprot,protl(*),ss_index(*)
      REAL*8  fstep,xp0(*),yp0(*),zp0(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER index(npm)
      REAL*8  gyr(npm)
      COMMON/rag1/ gyr,index

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,i1,j1,m,map,count,natom
      REAL*8  xpi,ypi,zpi,xc,yc,zc,rsp

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      map=0
      count=0
      DO i=1,nprot
         m=protl(map+1)
         IF(ss_index(protl(map+1+1)) .EQ. 1 .AND. m .NE. 0) THEN
            count=count+1
            gyr(count)=0.0D0
            index(count)=i
            natom=0
            DO j=1,m
               i1=protl(map+1+j)
               xpi=xp0(i1)
               ypi=yp0(i1)
               zpi=zp0(i1)
               DO k=j+1,m
                  j1=protl(map+1+k)
                  xc=xpi-xp0(j1)
                  yc=ypi-yp0(j1)
                  zc=zpi-zp0(j1)
                  rsp=DSQRT(xc*xc+yc*yc+zc*zc)
                  gyr(count)=gyr(count)+rsp
                  natom=natom+1
               END DO
            END DO
            gyr(count)=gyr(count)/DFLOAT(natom)
         END IF
         map=map+m+1
      END DO
      
      DO i=1,count,5
         m=4
         IF(count-i .LT. 4) m=count-i
         WRITE(kgyr,1000) 'T',fstep,(gyr(j),index(j),j=i,i+m)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1000  FORMAT(a1,f11.2,5(f12.5,1x,i5))
      RETURN
      END
