      SUBROUTINE search_clsth(protl,ma,nprot,concth,m1,mask1,mask2,nato
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <97/09/12 12:47:46 marchi>                             *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate the number of non-connected units in the system        *
*                                                                      *
*     protl  : List of the independent units                        O  *
*              >>INTEGER*4 protl(*)<<                                  *
*     mask1  : Total mask which indicates which atoms have been     W  *
*              gone through.                                           *
*              >>LOGICAL*4 mask1(*)<<                                  *
*     mask2  : Local mask reset to .TRUE. every loop through the    W  *
*              atoms of the units.                                     *
*              >>LOGICAL*4 mask1(*)<<                                  *
*     xp0    : Coordinates in simulation box coordinates          I/O  *
*     yp0      >>REAL*8  xp0(*),yp0(*),zp0(*)<<                        *
*     zp0                                                              *
*     nato   : Number of total atoms                                I  *
*              >>INTEGER*4 nato<<                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER m1
      INTEGER concth(m1,*),protl(*),ma,nprot,nato,iret
      LOGICAL mask1(*),mask2(*)
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ncount,m,i,nstart
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      ncount=1
      m=0
      DO i=1,nato
          mask1(i)=.TRUE.
      END DO
      ok=.TRUE.
      nprot=1
100   CONTINUE
      IF(m .NE. 0) THEN
         nprot=nprot+1
      END IF
      DO i=1,nato
          mask2(i)=.TRUE.
      END DO

*=======================================================================
*---  Peek the first atom never percolated upon and start new ----------
*---  percolation from there -------------------------------------------
*=======================================================================

      DO i=1,nato
         IF(mask1(i)) THEN
            nstart=i
            GOTO 200
         END IF
      END DO
      ok=.FALSE.
200   CONTINUE
      IF(ok) THEN

*=======================================================================
*---  Percolation start at nstart. On return mask2 ci .FALSE. for ------
*---  all atoms percolated upon ----------------------------------------
*=======================================================================

         CALL perc_bondh(nstart,concth,m1,mask2)

         m=0
         DO i=1,nato

*=======================================================================
*---  Accumulate atoms belonging to the cluster                 --------
*=======================================================================

           IF(.NOT. mask2(i)) THEN
               m=m+1
               protl(ncount+m)=i
               mask1(i)=.FALSE.
            END IF
         END DO

         IF(m .NE. 0) THEN
            IF(ncount .GE. ma) THEN
               iret=1
               errmsg=' From SEARCH_CLSTH: Physical dims.of CLSTHL'
     &              //' smaller than expected. Abort.'
               RETURN
            END IF 

*=======================================================================
*---  protl(ncount) contains the number of atoms of the cluster      ---
*=======================================================================

            protl(ncount)=m
            ncount=ncount+m+1
         END IF
         GOTO 100
      END IF
      nprot=nprot-1
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
