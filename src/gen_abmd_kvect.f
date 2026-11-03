      SUBROUTINE gen_abmd_kvect(vect,nvect,pvect)

************************************************************************
*   Time-stamp: <97/08/21 21:24:31 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Aug 21 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nvect
      REAL*8 vect,pvect(3,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  u1,u2,u3,w1,w2,w3,w,ranf
      external ranf

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nvect
100      CONTINUE
         u1=ranf()
         u2=ranf()
         u3=ranf()
         w1=1.0D0-2.0D0*u1
         w2=1.0D0-2.0D0*u2
         w3=1.0D0-2.0D0*u3
         w=w1**2+w2**2+w3**2
         IF(w .LT. 1.0D0) THEN
            pvect(1,i)=vect*w1/DSQRT(w)
            pvect(2,i)=vect*w2/DSQRT(w)
            pvect(3,i)=vect*w3/DSQRT(w)
         ELSE
            GOTO 100
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
