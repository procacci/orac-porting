      SUBROUTINE add1h(ax,length,x0,y0,z0)


************************************************************************
*                                                                      *
*     ADD1H will produce coordinates for 1 hydrogen relative           *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0, Y0, Z0 <<                               *
*                                                                      *
*---- Last update 07/08/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992               *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0,y0,z0,ax(3),length

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 norm,small
      INTEGER i
      DATA small/1.0D-4/

*==================== EXECUTABLE STATEMENTS ============================

      norm=0.0D0
      DO 10 i=1,3
          norm=norm+ax(i)**2
10    CONTINUE
      norm=DSQRT(norm)

      x0=ax(1)*length/norm+small
      y0=ax(2)*length/norm+small
      z0=ax(3)*length/norm+small

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
