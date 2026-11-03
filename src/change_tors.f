      SUBROUTINE change_tors(itor_ptype)

************************************************************************
*   Time-stamp: <97/02/08 11:13:58 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER itor_ptype

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parameters.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j

*----------------------- EXECUTABLE STATEMENTS ------------------------*


!=======================================================================
!------- Copy the three 1-dimensional parameter arrays to --------------
!---------------- a 3-dimensional array --------------------------------
!=======================================================================

      IF(itor_ptype .NE. 0) THEN
         DO i=1,lpito
            IF(itor_ptype .EQ. 1) THEN
               pito3(i)=1.0D0
            ELSE
               pito3(i)=-1.0D0
            END IF
         END DO
      END IF

      DO j=1,lpito
          pito(j,1)=pito1(j)
          pito(j,2)=pito2(j)
          pito(j,3)=pito3(j)
       END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
