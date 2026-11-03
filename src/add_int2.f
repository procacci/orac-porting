      SUBROUTINE add_int2(change,nato1,nato2,array1,array2,n1,n2,n0
     &     ,ndim1,nrep,iret)

************************************************************************
*   Time-stamp: <97/02/26 19:14:00 marchi>                             *
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

      INTEGER nato1,nato2,n0,n1,n2,iret,ndim1,nrep,offset
      INTEGER array1(n0,*),array2(n0,*)
      LOGICAL change

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,l,n,count

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
            DO n=1,n0
               IF(change) THEN
                  array1(n,l)=array2(n,i)+count
               ELSE 
                  array1(n,l)=array2(n,i)
               END IF
            END DO
         END DO
         count=count+nato2
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
