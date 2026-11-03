      SUBROUTINE get_sgr(line,cgroup)

************************************************************************
*   Time-stamp: <95/03/13 18:04:33 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 13 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      CHARACTER*80 line,cgroup      

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,count,nbeg(80),nend(80),length
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ok=.FALSE.
      count=0
      IF(line(1:1) .NE. ' ') THEN
         ok=.TRUE.
         nbeg(1)=1
         count=1
      END IF
      DO i=1,80
         cgroup(i:i)=' '
         IF(ok) THEN
            IF(line(i:i) .EQ. ' ') THEN
               nend(count)=i-1
               ok=.FALSE.
            END IF
         ELSE
            IF(line(i:i) .NE. ' ') THEN
               count=count+1
               nbeg(count)=i
               ok=.TRUE.
            END IF
         END IF
      END DO
      length=nend(count)-nbeg(4)+1
      cgroup(1:length)=line(nbeg(4):nend(count))

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
