      SUBROUTINE usort(unit,n1,label,ind,iret,errmsg)

************************************************************************
*                                                                      *
*     Add a new label to the list of residues.                         *
*                                                                      *
*     UNIT    :  List of labels for the residues and endgroups         *
*                >> character*1 UNIT(8,N1) <<                          *
*     N1      :  Physical column dimension of UNIT                     *
*     LABEL   :  Residue label to be compared with the ones in         *
*                the list.                                             *
*                >> character*1 LABEL(8) <<                            *
*     IND     :  Address of LABEL in UNIT.                             *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update 03/01/91 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley CA 1992                    *
*                                                                      *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nunit,ind,iret,n1
      CHARACTER*1 unit(8,*),label(8)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER m

*==================== EXECUTABLE STATEMENTS ============================

      ind=ind+1
      IF(ind .GT. n1) THEN
          iret=1
          errmsg=' In USORT: Number of residues exceed physical'//
     x           ' dimensions. Abort.'
          RETURN
      END IF
      DO m=1,8
          unit(m,ind)=label(m)
      END DO
      iret=0

*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
