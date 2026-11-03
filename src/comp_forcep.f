      SUBROUTINE comp_forcep(prt,st,oc,volume,pext)

************************************************************************
*   Time-stamp: <98/02/13 21:24:21 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Mar  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  prt(3,3),st(3,3),oc(3,3),volume,pext

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k
      REAL*8  sum,delta(3,3),aux(3,3)
      DATA delta/1.0D0,3*0.0D0,1.0D0,3*0.0D0,1.0D0/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO i=1,3
         DO j=1,3
            aux(i,j)=0.0D0
            DO k=1,3
               aux(i,j)=aux(i,j)+(prt(i,k)+st(i,k)-volume*pext*delta(i,k
     &              ))*oc(j,k)
            END DO
         END DO
      END DO
      DO i=1,3
         DO j=1,3
            prt(i,j)=aux(i,j)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
