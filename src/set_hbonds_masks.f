      SUBROUTINE set_hbonds_masks(nato,lacc,llacc,ldon,lldon,a_mask
     &     ,d_mask)

************************************************************************
*   Time-stamp: <97/07/10 18:55:55 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul 10 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato,llacc,lldon,lacc(2,*),ldon(2,*)
      LOGICAL a_mask(*),d_mask(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ia

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nato
         a_mask(i)=.FALSE.
         d_mask(i)=.FALSE.
      END DO
      DO i=1,llacc
         ia=lacc(1,i)
         a_mask(ia)=.TRUE.
      END DO
      DO i=1,lldon
         ia=ldon(1,i)
         d_mask(ia)=.TRUE.
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
