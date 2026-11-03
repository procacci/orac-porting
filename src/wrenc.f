      SUBROUTINE wrenc(kprint,line)

************************************************************************
*                                                                      *
*     Write a string of 78 characters between a star character.        *
*                                                                      *
*======================================================================*
*                                                                      *
*    Written by Massimo Marchi CEA/Centre d'Etudes Saclay FRANCE 1995  *
*                                                                      *
************************************************************************

*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint
      CHARACTER*80 line

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*1 char
      DATA char/'='/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(kprint,'(a1,a78,a1)') char,line(1:78),char

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
