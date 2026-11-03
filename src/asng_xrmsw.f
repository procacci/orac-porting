      SUBROUTINE asng_xrmsw(ss_point,ma,wca,whe,wbc,beta,mback,nbone)

************************************************************************
*   Time-stamp: <98/03/20 15:56:00 marchi>                             *
*                                                                      *
*   Assign weigths to atoms according to their type. CA's, Backbone    *
*   and heavy atoms are included.                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jun 29 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  wca(*),wbc(*),whe(*)
      INTEGER ma
      INTEGER nbone,mback(*),ss_point(ma,*)
      CHARACTER*7 beta(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,i,nmol,ntap

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nmol=ss_point(1,2)
      DO n=1,nmol
         i=ss_point(1+n,2)
         wca(i)=0.0D0
         wbc(i)=0.0D0
         whe(i)=0.0D0
      END DO
      ntap=ss_point(1,1)
      DO n=1,ntap
         i=ss_point(1+n,1)
         IF(beta(i) .EQ. 'ca') THEN
            wca(i)=1.0D0
         ELSE
            wca(i)=0.0D0
         END IF
      END DO
      DO i=1,nbone
         wbc(mback(i))=1.0D0
      END DO
      DO n=1,ntap
         i=ss_point(1+n,1)
         IF(beta(i)(1:1) .NE. 'h' .AND. beta(i)(1:2) .NE. 'cl' .AND.
     &        beta(i)(1:2) .NE. 'na') THEN
            whe(i)=1.0D0
         ELSE
            whe(i)=0.0D0
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
