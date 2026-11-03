      SUBROUTINE dupl_real1(array,na,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/03/14 15:24:00 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Mar 14 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER na,ndim,nrep,iret
      REAL*8  array(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=2,nrep
         DO m=1,na
            l=na*(n-1)+m
            IF(l .GT. ndim) THEN
               iret=1
               RETURN
            END IF
            array(l)=array(m)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
