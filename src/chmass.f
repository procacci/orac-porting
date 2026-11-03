      SUBROUTINE chmass(beta,mass,nato,hdmass)

************************************************************************
*                                                                      *
*     Change the mass of the hydrogens.                                *
*                                                                      *
*                                                                      *
*     MASS    :  List of atomic masses.               (INPUT/OUTPUT)   *
*                >> real*8   MASS(NTAP) <<                             *
*     NATO    :  Number of atoms making up the                (INPUT)  *
*                macromolecule.                                        *
*     HDMASS  :  New mass of the hydrogens                    (INPUT)  *
*                                                                      *
*                                                                      *
!---- Last update 02/16/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi  UC Berkeley 1991                      *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      INTEGER nato
      REAL*8  mass(*),hdmass
      CHARACTER*7 beta(*)

!-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j

!==================== EXECUTABLE STATEMENTS ============================

      DO i=1,nato
          IF(beta(i)(1:1) .EQ. 'h' .OR. beta(i)(1:1) .EQ. 'H') THEN
              mass(i)=hdmass
          END IF
      END DO

!=======================================================================

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
