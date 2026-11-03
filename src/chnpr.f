      SUBROUTINE chnpr(protl,ma,mb,nprot,concta,m1,mask1,mask2,nato
     &     ,iret,errmsg)

************************************************************************
*                                                                      *
*     Calculate the number of indipendent units which compose          *
*     the system.                                                      *
*                                                                      *
*     protl  : List of the independent units                        O  *
*              >>INTEGER*4 protl(*)<<                                  *
*     mask1  : Total mask which indicates which atoms have been     W  *
*              gone through.                                           *
*              >>LOGICAL* mask1(*)<<                                  *
*     mask2  : Local mask reset to .TRUE. every loop through the    W  *
*              atoms of the units.                                     *
*              >>LOGICAL*4 mask1(*)<<                                  *
*     xp0    : Coordinates in simulation box coordinates          I/O  *
*     yp0      >>REAL*8  xp0(*),yp0(*),zp0(*)<<                        *
*     zp0                                                              *
*     nato   : Number of total atoms                                I  *
*              >>INTEGER*4 nato<<                                      *
*                                                                      *
*---- Last update 27/10/93 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CEA, Saclay 1993                       *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      IMPLICIT none
 
      INTEGER m1
      INTEGER protl(*),ma,mb,nprot,nato,iret,concta(m1,*)
      LOGICAL mask1(*),mask2(*)
      CHARACTER*80 errmsg

      INTEGER ncount,m,i,nstart
      LOGICAL ok
 
      ncount=1
      DO i=1,nato
          mask1(i)=.TRUE.
      END DO
      ok=.TRUE.
      nprot=0
100   CONTINUE
      nprot=nprot+1
      DO i=1,nato
          mask2(i)=.TRUE.
      END DO
      DO i=1,nato
          IF(mask1(i)) THEN
              nstart=i
              GOTO 200
          END IF
      END DO
      ok=.FALSE.
200   CONTINUE
      IF(ok) THEN
          CALL recpr(nstart,concta,m1,mask2)
          m=0
          DO i=1,nato
              IF(.NOT. mask2(i)) THEN
                  m=m+1
                  protl(ncount+m)=i
                  mask1(i)=.FALSE.
              END IF
          END DO
          IF(m .NE. 0) THEN
             IF(ncount .GE. ma) THEN
                iret=1
                errmsg=' In CHNPR: Physical dimensions of protl smaller'
     &               //' than expected. Abort.'
                RETURN
             END IF 
             protl(ncount)=m
             ncount=ncount+m+1
          END IF
          GOTO 100
      END IF
      nprot=nprot-1
    
      RETURN
      END
