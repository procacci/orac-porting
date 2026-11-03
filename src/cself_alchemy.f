      SUBROUTINE cself_alchemy(ss_index,nato,alphal,kcut,charge,charge0
     &     ,self_slt,dwrk)

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
      REAL*8  charge(*),charge0(*),alphal,self_slt,self_slv,kcut,dwrk

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,n
      REAL*8  fsrtal(2)
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,alphar,en0
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the self term --------------------------------------------
*=======================================================================

      fsrtal(1)=0.0D0
      en0=0.d0
      DO n=1,nato
         i=ss_index(n)
         if(i.eq.1) THEN 
           fsrtal(i)=fsrtal(i)+charge(n)**2
           en0=en0+charge0(n)**2
         END IF
      END DO
      alphar=kcut/(2.0D0*alphal) 
      qt=1.0D0/(1.0D0+qp*alphar)
      expcst=exp(-alphar*alphar)
      erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
      fsrtal(1)=-fsrtal(1)*alphal/DSQRT(pi)*(1.0d0 - erfcst)
      en0=-en0*alphal/DSQRT(pi)*(1.0d0 - erfcst)
      self_slt=fsrtal(1)
      dwrk=fsrtal(1)-en0

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
