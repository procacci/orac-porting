      SUBROUTINE add_concta(nato1,nato2,array1,array2,n1,n2,n0
     &     ,ndim1,ndim2,nrep,iret)

************************************************************************
*   Time-stamp: <97/02/11 17:41:40 marchi>                             *
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

      INTEGER nato1,nato2,n0,n1,n2,iret,ndim1,ndim2,nrep,offset
      INTEGER array1(ndim1,*),array2(ndim2,*)
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
            array1(l,1)=array2(i,1)
            DO n=2,n0
               array1(l,n)=array2(i,n)+count
            END DO
         END DO
         count=count+nato2
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
