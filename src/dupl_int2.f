      SUBROUTINE dupl_int2(change,offset,narray,na,n1,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/11/02 14:30:57 marchi>                             *
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
      
      INTEGER ndim,n1,narray(n1,*),na,nrep,iret,offset
      LOGICAL change

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,k,count

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
            DO k=1,n1
               IF(change) THEN
                  narray(k,l)=narray(k,m)+count
               ELSE
                  narray(k,l)=narray(k,m)
               END IF
            END DO
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
