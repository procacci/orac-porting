      SUBROUTINE dupl_char1(carray,na,nb,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/03/22 15:16:52 marchi>                             *
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
      
      INTEGER na,nb,ndim,nrep,iret
      CHARACTER*1 carray(nb,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,k

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=2,nrep
         DO m=1,na
            l=na*(n-1)+m
            IF(l .GT. ndim) THEN
               iret=1
               RETURN
            END IF
            DO k=1,nb
               carray(k,l)=carray(k,m)
            END DO
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
