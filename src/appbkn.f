      SUBROUTINE appbkn(label,natom,cback,nback,ibackp,shifta,shiftb,
     x           iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine appends a list of atom labels, types,            *
*     charges, residue pointers to existing lists.                     *
*                                                                      *
*                                                                      *
*     LABEL     :  List of labels for each atom of the residue         *
*                  >> character*7 LABEL(NATOM) <<             (INPUT)  *
*     NATOM     :  Number of atoms of the current residue     (INPUT)  *
*                  >> integer NATOM <<                                 *
*     CBACK     :  List of backbone atoms                     (INPUT)  *
*                  >> character*7 CBACK(NBACK) <<                      *
*     NBACK     :  Number of backbone atoms                   (INPUT)  *
*                  >> integer NBACK <<                                 *
*     IBACKP    :  List of atoms of the backbone             (OUTPUT)  *
*                  >> integer IBACKP(*) <<                             *
*     SHIFTA    :  Number of atoms currently appended         (INPUT)  *
*                  >> integer NSHIFT <<                                *
*     SHIFTB    :  Number of backbone atoms currently         (INPUT)  *
*                  appended                                            *
*                  >> integer NSHIFT <<                                *
*     IRET      :  Return code                               (OUTPUT)  *
*                  >> integer IRET <<                                  *
*                  IRET=0   : Normal termination                       *
*                  IRET=1   : Fatal return code                        *
*     ERRMSG    :  Error message                             (OUTPUT)  *
*                  >> character*80 errmsg <<                           *
*                                                                      *
*---- Last update 03/31/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER natom,nback,shifta,shiftb
      INTEGER ibackp(*),iret
      CHARACTER*7 cback(*),label(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m
      LOGICAL err,ok

*==================== EXECUTABLE STATEMENTS ============================


*-----------------------------------------------------------------------
*------------- Append backbone atoms to the list -----------------------
*-----------------------------------------------------------------------

      n=1
      err=.false.
10    IF((.NOT.err).AND.(n.le.nback)) THEN
          m=1
          ok=.false.
20        IF((.NOT.ok).AND.(m.le.natom)) THEN
              IF(cback(n).EQ.label(m)) THEN
                  ibackp(n+shiftb)=m+shifta
                  ok=.true.
              END IF
              m=m+1
              GOTO 20
          END IF
          n=n+1
          err=.NOT.ok
          GOTO 10
      END IF


*-------------- Check if an error has occurred ------------------------

      IF(err) THEN
          iret=1
          errmsg=' Fatal error in APPBKN : backbone atom NOT FOUND ** '
      ELSE
          iret=0
          errmsg='  '
      END IF
*-----------------------------------------------------------------------


      RETURN
      END
