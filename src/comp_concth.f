      SUBROUTINE comp_concth(ntap,beta,concta,concth,m1)

************************************************************************
*   Time-stamp: <97/03/28 12:01:09 marchi>                             *
*                                                                      *
*                                                                      *
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

      INTEGER m1,ntap
      INTEGER concta(m1,*),concth(m1,*)
      CHARACTER*7 beta(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER m,n,i,j,k

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,ntap
         m=concta(i,1)
         concth(i,1)=0
         n=0
         DO j=1,m
            k=concta(i,j+1)
            IF(beta(k)(1:1) .EQ. 'h' .OR. beta(i)(1:1) .EQ. 'h ') THEN
               n=n+1
               concth(i,n+1)=k
            END IF
         END DO
         concth(i,1)=n
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
