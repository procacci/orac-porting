      SUBROUTINE check_topology(label,list,number,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/08 20:12:23 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      CHARACTER*1 label
      INTEGER number,iret,list(*)
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,i1
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ok=.TRUE.
      DO i1=1,list(1)
         i=list(1+i1)
         IF(number .LT. i) ok=.FALSE.
         IF(i .LT. 0) ok=.FALSE.
      END DO
      IF(.NOT. ok) THEN
         iret=1
         IF(label .EQ. 'B') THEN
            errmsg=
     &' While writing bond lengths. One or more bonds not in the list. '
         ELSE IF(label .EQ. 'D') THEN
            errmsg=
     &' While writing bond angles. One or more bendings not '
     &//'in the list. '
         ELSE IF(label .EQ. 'P') THEN
            errmsg=
     &' While writing proper torsions. One or more torsions '
     &//'not in the list. '
         ELSE IF(label .EQ. 'I') THEN
            errmsg=
     &' While writing improper torsions. One or more torsions '
     &//'not in the list. '         
         END IF
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
