      SUBROUTINE find_thermos(nato,cnstpp,nato_slt,nato_slv,nmol_slt
     &     ,nmol_slv,cnstpp_slv,nprot,cpress,isostress,ndf_thermos)

************************************************************************
*   Time-stamp: <98/03/18 08:19:56 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Apr  5 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nato,cnstpp,nmol_slv,nmol_slt,nato_slv,nato_slt
     &     ,cnstpp_slv,nprot,ndf_thermos(*)
      LOGICAL cpress,isostress

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nfree_slv,nfree_slt,cnstpp_slt
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ndf_thermos(1)=3*nprot

      cnstpp_slt=cnstpp-cnstpp_slv
      IF(cpress) THEN
         IF(isostress) THEN
            ndf_thermos(1)=ndf_thermos(1)+1
         ELSE
            ndf_thermos(1)=ndf_thermos(1)+6
         END IF
      END IF

      nfree_slv=nato_slv*3-3*nmol_slv-cnstpp_slv
      nfree_slt=nato_slt*3-3*nmol_slt-cnstpp_slt
      ndf_thermos(2)=nfree_slt
      ndf_thermos(3)=nfree_slv

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
