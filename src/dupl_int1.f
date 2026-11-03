      SUBROUTINE dupl_int1(change,offset,narray,na,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/03/22 17:34:33 marchi>                             *
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
      
      INTEGER narray(*),na,ndim,nrep,iret,offset
      LOGICAL change

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,count

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO n=2,nrep
         count=count+offset
         DO m=1,na
            l=na*(n-1)+m
            IF(l .GT. ndim) THEN
               iret=1
               RETURN
            END IF
            IF(change) THEN
               narray(l)=narray(m)+count
            ELSE
               narray(l)=narray(m)
            END IF
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
