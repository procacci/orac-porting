      SUBROUTINE comp_rmsq(time,rms_disp,xpc,ypc,zpc,tot_rms_disp
     &     ,dsq,tnormb)

************************************************************************
*   Time-stamp: <97/12/08 15:41:06 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Nov 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER time
      REAL*8  xpc(*),ypc(*),zpc(*),rms_disp(*),tot_rms_disp(*)
     &     ,dsq(0:*),tnormb(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,no2,j
      REAL*8  sumsq,fcorr,drx,dry,drz

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,time
         tnormb(i) = 1.0D0/DFLOAT(time-i+1)
      END DO
      sumsq=0.0D0
      DO i=1,time
         dsq(i)=xpc(i)**2+ypc(i)**2+zpc(i)**2
         sumsq=sumsq+2.0D0*dsq(i)
      END DO
      dsq(0)=0.0D0
      dsq(time+1)=0.0D0
      DO i=1,time-1
         sumsq=sumsq-dsq(i-1)-dsq(time-i+2)
         tot_rms_disp(i) = tot_rms_disp(i)+tnormb(i)*(sumsq - 2.0D0
     &        *rms_disp(i))
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
