      SUBROUTINE fpbond_adjust_bonds(lbnd,lbond,xp0,yp0,zp0,pota,
     &     potb,ubond)

************************************************************************
*   Time-stamp: <97/02/28 12:10:25 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 23 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER lbond,lbnd(2,*)
      REAL*8  ubond
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb
      REAL*8  xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21

*==================== EXECUTABLE STATEMENTS ============================

      ubond=0.0D0

*----- bonded interactions: bond stretching

      DO i=1,lbond
          la=lbnd(1,i)
          lb=lbnd(2,i)

          xr1=xp0(la)
          yr1=yp0(la)
          zr1=zp0(la)
          xr2=xp0(lb)
          yr2=yp0(lb)
          zr2=zp0(lb)

          x21=xr2-xr1
          y21=yr2-yr1
          z21=zr2-zr1
          rs21=DSQRT(x21**2+y21**2+z21**2)
          ubond=ubond+potb(i)*(rs21-pota(i))**2
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
