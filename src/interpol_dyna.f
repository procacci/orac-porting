      SUBROUTINE interpol_dyna(a,b,ccoeff,cder,fc,nstep,dyna)

************************************************************************
*   Time-stamp: <98/03/04 23:05:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Mar  4 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstep
      REAL*8  a,b,ccoeff(*),cder(*),fc(*),dyna,chebev
      EXTERNAL chebev
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER l

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL chebft_dyn(ccoeff,nstep,fc)
      CALL chder_dyn(a,b,ccoeff,cder,nstep)
      dyna=chebev(a,b,cder,nstep,0.0D0)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
