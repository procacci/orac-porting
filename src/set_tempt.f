      SUBROUTINE set_tempt(n,qmass,etap,temph,t)

************************************************************************
*   Time-stamp: <97/04/06 12:20:04 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Apr  6 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n
      REAL*8  qmass(*),etap(*),temph,t

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  ranf,tto,sig,tvel,u1,u2,t1
      EXTERNAL ranf

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO i=1,n
         IF(DABS(qmass(i)) .GT. 1.0D-7) THEN
            tvel=boltz*t/(unite*qmass(i))
            sig=DSQRT(tvel)
            u1=ranf()
            u2=ranf()
            t1=DSQRT(-2.0d0*DLOG(u1))*DCOS(2.0d0*pi*u2)
            etap(i)=t1*sig
         ELSE
            etap(i)=0.0D0
         END IF
      END DO
      tto=0.0D0
      DO i=1,n
         tto=tto+0.50D0*qmass(i)*tto/gascon
      END DO
      tto=tto*efact
      temph=(2.0D0/3.0D0)*tto/gascon

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
