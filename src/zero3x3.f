      SUBROUTINE zero3x3(mat)

************************************************************************
*                                                                      *
*     ZERO will zero averages for the run.                             *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  mat(3,3),zero
      PARAMETER (zero=0.0D0)

*==================== EXECUTABLE STATEMENTS ============================

      mat(1,1)=zero
      mat(1,2)=zero
      mat(1,3)=zero
      mat(2,1)=zero
      mat(2,2)=zero
      mat(2,3)=zero
      mat(3,1)=zero
      mat(3,2)=zero
      mat(3,3)=zero

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
