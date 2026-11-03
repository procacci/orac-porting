      SUBROUTINE dupl_concta(offset,narray,na,n1,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/11/02 14:28:58 marchi>                             *
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
      
      INTEGER ndim,narray(ndim,*),na,n1,nrep,iret,offset

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
            narray(l,1)=narray(m,1)
            DO k=2,n1
               narray(l,k)=narray(m,k)+count
            END DO
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
