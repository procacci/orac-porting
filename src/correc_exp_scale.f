      SUBROUTINE correc_exp_scale(cpress,slt_exist,slv_exist,ma,mb,etap
     &     ,tm4,nprot,vcax,vcay,vcaz,vpx,vpy,vpz,vco)

************************************************************************
*   Time-stamp: <98/02/12 00:05:58 marchi>                             *
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
      
      INTEGER ma(*),mb(*),nprot
      LOGICAL cpress,slv_exist,slt_exist
      REAL*8  tm4,vcax(*),vcay(*),vcaz(*),vpx(*),vpy(*),vpz(*),vco(3,3)
     &     ,etap(*)

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

c----------  Correct cofm and barostat velocities using thermostat 1 ---
c --
      CALL correc_scale(1,nprot,vcax,vcay,vcaz,etap(1),tm4)
               
      IF(cpress) THEN
         CALL correc_scale_co(9,vco,etap(1),tm4)
      END IF
      
c----------  Correct solute velocities using thermostat 2 --------------
c --

      
      IF(slt_exist) CALL correc_scale_2(ma,vpx,vpy,vpz,etap(2),tm4)
      
c----------  Correct solvent velocities using thermostat 3 -------------
c --
      IF(slv_exist) CALL correc_scale_2(mb,vpx,vpy,vpz,etap(3),tm4)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
