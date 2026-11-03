      SUBROUTINE add_int1(change,nato1,nato2,n1,n2,array1,array2,ndim1
     &     ,nrep,iret)

************************************************************************
*   Time-stamp: <97/02/11 09:25:40 marchi>                             *
*                                                                      *
*   Add content of array2 to array1                                    *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Feb 10 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato1,nato2,n1,n2,iret,ndim1,nrep,offset
      INTEGER array1(*),array2(*)
      LOGICAL change

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,l,count

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(ndim1 .LT. n1+n2*nrep) THEN
         iret=1
         RETURN
      END IF
      
      count=nato1
      l=n1
      DO j=1,nrep
         DO i=1,n2
            l=l+1
            IF(change) THEN
               array1(l)=array2(i)+count
            ELSE 
               array1(l)=array2(i)
            END IF
         END DO
         count=count+nato2
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
