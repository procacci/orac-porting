      SUBROUTINE check_read_columns(nbuffer,length_run,atom_record
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <97/11/26 15:00:23 marchi>                             *
*                                                                      *
*   Determine how many columns can be read at the time given           *
*   a dimension of nbuffer. nblock is the number of atoms that         *
*   can be read at the time.                                           *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Nov 25 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nbuffer,length_run,atom_record,iret
      CHARACTER*80 errmsg

*---------------------------LOCAL VARIABLES ---------------------------*

      CHARACTER*8 string
      INTEGER number,time

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      time=nbuffer/atom_record
      number=atom_record*length_run
      WRITE(string,'(i8)') number
      IF(time .LT. length_run) THEN
         iret=1
         errmsg=
     &        'Buffer length not sufficient: increase to at least:'
     &        / /string
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
