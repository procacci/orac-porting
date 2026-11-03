      SUBROUTINE collision(ntap,vpx0,vpy0,vpz0,massb,nutime,t,time)

************************************************************************
*   Time-stamp: <98/01/31 19:26:39 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Sep 13 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER ntap
      REAL*8  vpx0(*),vpy0(*),vpz0(*),massb(*),nutime,t,time

*------------------------- LOCAL VARIABLES ----------------------------*

      LOGICAL near0
      REAL*8  tvel,u1,u2,u3,u4,ranf,exp_dev,vpx,vpy,vpz,t1,t2,t3,sig
      INTEGER time_colp(m1), i
      DATA time_colp/m1*0/
      SAVE time_colp
      external ranf

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,ntap
         IF(time_colp(i) .EQ. 0) THEN
200         exp_dev=ranf()
            IF(exp_dev .EQ. 0.0D0) GOTO 200
            exp_dev=-nutime*DLOG(exp_dev)
            time_colp(i)=IDINT(exp_dev/time)
            tvel=0.0D0
            IF(.NOT. near0(massb(i))) tvel=boltz*t/(unite*massb(i))
            sig=DSQRT(tvel)
            u1=ranf()
            u2=ranf()
            u3=ranf()
            u4=ranf()
            IF(u1 .EQ. 0.0D0) THEN
               u1=1.0D-15
            END IF
            IF(u2 .EQ. 0.0D0) THEN
               u2=1.0D-15
            END IF
            IF(u3 .EQ. 0.0D0) THEN
               u3=1.0D-15
            END IF
            IF(u4 .EQ. 0.0D0) THEN
               u4=1.0D-15
            END IF
            t1=DSQRT(-2.0d0*DLOG(u1))*DCOS(2.0d0*pi*u2)
            t2=DSQRT(-2.0d0*DLOG(u1))*DSIN(2.0d0*pi*u2)
            t3=DSQRT(-2.0d0*DLOG(u3))*DCOS(2.0d0*pi*u4)
            vpx0(i)=t1*sig
            vpy0(i)=t2*sig
            vpz0(i)=t3*sig
         ELSE
            time_colp(i)=time_colp(i)-1
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
