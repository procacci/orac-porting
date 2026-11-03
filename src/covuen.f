      SUBROUTINE covuen(alfa,emin,emax)

************************************************************************
*                                                                      *
*     Convert folding parameters to program units.                     *
*                                                                      *
*---- Last update 04/01/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley 1991                       *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  alfa,emin,emax

!==================== EXECUTABLE STATEMENTS ============================

      alfa=alfa*unite*avogad/(1000.0D0)
      emin=emin*1000.0D0/(unite*avogad)
      emax=emax*1000.0D0/(unite*avogad)

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
