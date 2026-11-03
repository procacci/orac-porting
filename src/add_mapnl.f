      SUBROUTINE add_mapnl(nato1,nato2,n1,n2,mapnl,mapnl_slv,ndim1
     &     ,nrep,iret)

************************************************************************
*   Time-stamp: <97/02/11 09:44:25 marchi>                             *
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

      INTEGER nato1,nato2,n1,n2,iret,ndim1,nrep
      INTEGER mapnl(*),mapnl_slv(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER j,l,m,count,map1,map2,ma

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      map2=0
      DO m=1,nato1
         ma=mapnl(map2+1)
         map2=map2+ma+1
      END DO
      map1=0
      DO m=1,nato2
         ma=mapnl_slv(map1+1)
         map1=map1+ma+1
      END DO

      IF(ndim1 .LT. map2+map1*nrep) THEN
         iret=1
         RETURN
      END IF

      count=nato1
      DO j=1,nrep
         map1=0
         DO m=1,nato2
            ma=mapnl_slv(map1+1)
            mapnl(map2+1)=ma
            DO l=1,ma
               mapnl(map2+1+l)=mapnl_slv(map1+1+l)+count
            END DO
            map1=map1+ma+1
            map2=map2+ma+1
         END DO
         count=count+nato2
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
