      SUBROUTINE dupl_mapnl(offset,mapnl,na,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/03/24 17:58:17 marchi>                             *
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
      
      INTEGER mapnl(*),na,ndim,nrep,iret,offset

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,count,map1,map2,ma

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      map2=0
      DO m=1,na
         ma=mapnl(map2+1)
         map2=map2+ma+1
      END DO
      
      count=0
      DO n=2,nrep
         count=count+offset
         map1=0
         DO m=1,na
            ma=mapnl(map1+1)
            mapnl(map2+1)=ma
            DO l=1,ma
               mapnl(map2+1+l)=mapnl(map1+1+l)+count
            END DO
            map1=map1+ma+1
            map2=map2+ma+1
            IF(map2 .GT. ndim) THEN
               iret=1
               RETURN
            END IF
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
