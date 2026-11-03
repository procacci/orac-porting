      SUBROUTINE min_pack(nato,map,nmax)

************************************************************************
*   Time-stamp: <95/03/27 18:57:31 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 27 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato,nmax,map(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,count

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      count=0
      DO n=1,nato
         m=map(count+1)
         count=count+m+1
      END DO
      nmax=count

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
