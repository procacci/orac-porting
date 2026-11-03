      SUBROUTINE csort(unit,nunit,label,ind,iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine compares the label of a residue or endgroup      *
*     with a list of labels and finds its address in that list.        *
*                                                                      *
*     UNIT    :  List of labels for the residues and endgroups         *
*                >> character*1 UNIT(8,N1) <<                          *
*     NUNIT   :  Number of labels in the list.                         *
*     N1      :  Physical column dimension of UNIT                     *
*     LABEL   :  Residue label to be compared with the ones in         *
*                the list.                                             *
*                >> character*1 LABEL(8) <<                            *
*     IND     :  Address of LABEL in UNIT.                             *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update 09/04/89 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nunit,ind,iret
      CHARACTER*1 unit(8,*),label(8)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m,map

*==================== EXECUTABLE STATEMENTS ============================

      n=0
      DO 10 map=1,nunit
          n=0
          DO 20 m=1,8
              IF(label(m).EQ.unit(m,map)) THEN
                  n=n+1
              END IF
20        CONTINUE
          IF(n.EQ.8) GOTO 100
10    CONTINUE
      iret=1
      errmsg=' Could not find RESIDUE ' // label(1)// label(2)//
     &     label(3) //label(4) // label(5)// label(6) 
     &     //label(7)// label(8) 
      RETURN

100   CONTINUE
      iret=0
      ind=map

*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
