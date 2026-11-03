      SUBROUTINE get_velocities(x,y,vx,length,fstep)

************************************************************************
*   Time-stamp: <97/12/04 15:11:59 marchi>                             *
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

      REAL*8  fstep,x(*),y(4,*),vx(*)
      INTEGER length

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,divide,bin
      REAL*8  time,dt,h

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      dt=x(2)-x(1)
      DO i=1,length-2
         time=fstep*(i-1)
         bin=time/dt+1
         h=time-x(bin)
         vx(i)=y(2,bin)+h*(y(3,bin)+h*y(4,bin)/3.0D0)/2.0D0
      END DO
      vx(length-1)=0.0D0
      vx(length)=0.0D0

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      
      RETURN
      END
