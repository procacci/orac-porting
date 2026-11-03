      FUNCTION f1dim_adjust_bonds(x,xp0,yp0,zp0,fpx,fpy
     &     ,fpz)

************************************************************************
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),x

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i
      REAL*8  ubond,f1dim_adjust_bonds

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------
  

*==================== EXECUTABLE STATEMENTS ============================


      DO i=1,ntap
         xp0(i)=xp0(i)+x*fpx(i)
         yp0(i)=yp0(i)+x*fpy(i)
         zp0(i)=zp0(i)+x*fpz(i)
      END DO
      CALL fpbond_adjust_bonds(lcnstr,lconstr,xp0,yp0,zp0,potbo_cnst(1,2
     &     ),potbo_cnst(1,1),ubond)
      
      DO i=1,ntap
         xp0(i)=xp0(i)-x*fpx(i)
         yp0(i)=yp0(i)-x*fpy(i)
         zp0(i)=zp0(i)-x*fpz(i)
      END DO
      f1dim_adjust_bonds=ubond

      RETURN
      END
