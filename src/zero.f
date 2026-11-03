      SUBROUTINE zero(avp,avs,m1,m2)

************************************************************************
*                                                                      *
*     ZERO will zero averages for the run.                             *
*                                                                      *
*     AVP     :  Solute averages.                                 (I)  *
*                >> real*8 AVP(M1) <<                                  *
*     AVS     :  Solvent averages.                                (I)  *
*                >> real*8 AVS(M1) <<                                  *
*     M1      :  Physical dimension of AVP.                       (I)  *
*     M2      :  Physical dimension of AVS.                       (I)  *
*                                                                      *
*     XSM     :  Displacement from the initial position           (I)  *
*     YSM        of each solvent molecule.                             *
*     ZSM        >> real*8 XSM(L2), ... <<                             *
*     L2      :  Physical dimension of XSM/YSM/ZSM.               (I)  *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS ZEROA.                                                 *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER m1,m2
      REAL*8  avp(m1),avs(m2)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,n

*==================== EXECUTABLE STATEMENTS ============================

      DO 30 i=1,m1
          avp(i)=0.0D+0
30    CONTINUE
      DO 40 i=1,m2
          avs(i)=0.0D+0
40    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
