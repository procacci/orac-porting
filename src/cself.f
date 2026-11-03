      SUBROUTINE cself(ss_index,nato,alphal,kcut,charge,self_slt
     &     ,self_slv)

************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald self term.                                     *
*                                                                      *
*---- Last update 06/12/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS  NONE                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================
      
      use unit
      
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ss_index(*),nato
      REAL*8  charge(*),alphal,self_slt,self_slv,kcut

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,n
      REAL*8  fsrtal(2)
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,alphar
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the self term --------------------------------------------
*=======================================================================

      fsrtal(1)=0.0D0
      fsrtal(2)=0.0D0
      DO n=1,nato
         i=ss_index(n)
         fsrtal(i)=fsrtal(i)+charge(n)**2
      END DO
      DO n=1,2
         alphar=kcut/(2.0D0*alphal)
         qt=1.0D0/(1.0D0+qp*alphar)
         expcst=exp(-alphar*alphar)
         erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
         fsrtal(n)=-fsrtal(n)*alphal/DSQRT(pi)*(1.0d0 - erfcst)
      END DO
      self_slt=fsrtal(1)
      self_slv=fsrtal(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
