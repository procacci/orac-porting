      SUBROUTINE find_length_run(stop_anl)

************************************************************************
*   Time-stamp: <97/12/07 14:55:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Nov 26 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER stop_anl

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER rec,uniit
      LOGICAL ok
      REAL*4 faux
      LOGICAL near0

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*----- Read box and timestep -------------------------------------------
*=======================================================================

      uniit=kwrite_dump(1)
      ok=.TRUE.
      rec=0
100   CONTINUE
         rec=rec+1
         READ(unit=uniit,rec=rec,ERR=200) faux
         IF(ABS(faux) .LT. 1.0e-6) GOTO 200
         GOTO 100
200   CONTINUE 
      stop_anl=rec-1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
