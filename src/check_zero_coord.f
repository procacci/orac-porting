      SUBROUTINE check_zero_coord(xp,yp,zp,nato,iret,errmsg)

************************************************************************
*   Time-stamp: <98/03/18 16:21:17 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Mar 18 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 xp(*),yp(*),zp(*)
      INTEGER nato,iret
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,count
      REAL*8  aux

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      count=0
      DO i=1,nato
         aux=xp(i)**2+yp(i)**2+zp(i)**2
         IF(DABS(aux) .LT. 1.0D-10) count=count+1
      END DO

      IF(count .GT. 1) THEN
         errmsg=
     &   ' WARNING: at least two coordinates in trajectory are ' / /
     &   'identically zero.'
         iret=1
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
