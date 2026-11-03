      SUBROUTINE comp_gmgp(dt,co,oc,vco,gmgp)

************************************************************************
*   Time-stamp: <97/04/08 10:50:21 marchi>                             *
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

      REAL*8  co(3,3),oc(3,3),vco(3,3),gmgp(3,3),dt

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  one(3,3),two(3,3),wk(9),sum,eigvl(3)
     &     ,eigvc(3,3)
      INTEGER ier,i,j,k

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO i=1,3
         DO j=1,3
            sum=0.0D0
            DO k=1,3
               sum=sum+vco(k,i)*co(k,j)+vco(k,j)*co(k,i)
            END DO
            one(i,j)=sum
         END DO
      END DO

      DO i=1,3
         DO j=1,3
            sum=0.0D0
            DO k=1,3
               sum=sum+oc(i,k)*oc(j,k)
            END DO
            two(i,j)=sum
         END DO
      END DO
      
      DO i=1,3
         DO j=1,3
            sum=0.0D0
            DO k=1,3
               sum=sum+one(i,k)*two(k,j)
            END DO
            gmgp(i,j)=sum
         END DO
      END DO

      CALL eigrs(gmgp,3,12,eigvl,eigvc,3,wk,ier)

      DO i=1,3
         DO j=1,3
            sum=0.0D0
            DO k=1,3
               sum=sum+eigvc(i,k)*eigvc(j,k)*DEXP(-dt*eigvl(k))
            END DO
            gmgp(i,j)=sum
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
