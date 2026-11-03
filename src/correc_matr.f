      SUBROUTINE correc_matr(tm2,co,oc,vco,gmgp,vcax,vcay,vcaz,nprot)

************************************************************************
*   Time-stamp: <98/02/10 19:20:49 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprot
      REAL*8 tm2,gmgp(3,3),co(3,3),oc(3,3),vco(3,3),vcax(*),vcay(*)
     &     ,vcaz(*)

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*


c--- Compute $e^{G^{-1} \dot G}$ ---------------------------------------
               
      CALL comp_gmgp(tm2,co,oc,vco,gmgp)
      
c--- Correct velocities of the cofm by gmgp ----------------------------
      
      CALL correc_gmgp(vcax,vcay,vcaz,gmgp,nprot)
               
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
