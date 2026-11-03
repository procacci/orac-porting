      SUBROUTINE openf(kunit,file,form1,status1,recl)

************************************************************************
*   Time-stamp: <98/03/09 16:39:49 marchi>                             *
*                                                                      *
*   Connect a file to a unit                                           *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Jul 17 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      CHARACTER*1  form1(*),status1(*)
      CHARACTER*80 file
      INTEGER      kunit,recl

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER unit,i
      CHARACTER*80 form,status
      SAVE unit
      DATA unit/10/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,80
         form(i:i)=' '
         status(i:i)=' '
      END DO
      IF(form1(1) .EQ. 'F') form='FORMATTED'
      IF(form1(1) .EQ. 'U') form='UNFORMATTED'
      IF(status1(1) .EQ. 'U') status='UNKNOWN'
      IF(status1(1) .EQ. 'N') status='NEW'
      IF(status1(1) .EQ. 'O') status='OLD'
      kunit=unit
      IF(recl .EQ. 0) THEN
         OPEN(unit=unit,file=file,form=form,status=status)
      ELSE
         OPEN(unit=unit,file=file,access='DIRECT',form=form,status
     &        =status,recl=recl)
      END IF
      IF(form1(1) .EQ. 'F') THEN
         REWIND unit
      END IF
      unit=unit+1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
