      SUBROUTINE setup_skin_shell(coupl_grp,coupl_mol,h,l,m,n1,n0)

************************************************************************
*   Time-stamp: <98/02/11 09:36:31 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL coupl_grp,coupl_mol,h,l,m,n1,n0

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(coupl_mol) THEN
         h=.FALSE.
         l=.FALSE.
         m=.TRUE.
         n1=.FALSE.
         n0=.FALSE.
      ELSE IF(coupl_grp) THEN
         h=.FALSE.
         l=.FALSE.
         m=.FALSE.
         n1=.FALSE.
         n0=.TRUE.
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
