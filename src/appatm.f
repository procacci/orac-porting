      SUBROUTINE appatm(label,type,charg1,mass1,nres,atm,natom,alpha
     &     ,beta,qchar1,qmass1,nbeta,ntd,nshift)

************************************************************************
*                                                                      *
*     This subroutine appends a list of atom labels, types,            *
*     charges and residue indexes to existing lists.                   *
*                                                                      *
*                                                                      *
*     LABEL     :  List of labels for each atom of the residue         *
*                  >> character*7 LABEL(NATOM) <<             (INPUT)  *
*     TYPE      :  List of label types for each atom          (INPUT)  *
*                  >> character*7 TYPE(NATOM) <<                       *
*     CHARG1    :  List of charges for each atom              (INPUT)  *
*                  >> real*8 CHARG1(NATOM) <<                          *
*     NRES      :  Type of the current residue                (INPUT)  *
*                  >> integer NRES <<                                  *
*     ATM       :  Number of the current residue              (INPUT)  *
*                  >> integer ATM <<                                   *
*     NATOM     :  Number of atoms of the current residue     (INPUT)  *
*                  >> integer NATOM <<                                 *
*     ALPHA     :  List of labels to which LABEL             (OUTPUT)  *
*                  has to be appended                                  *
*                  >> character*7 ALPHA(*) <<                          *
*     BETA      :  List of labels to which TYPE              (OUTPUT)  *
*                  has to be appended                                  *
*                  >> character*7 BETA(*) <<                           *
*     QCHAR1    :  List of charges to which CHARG1           (OUTPUT)  *
*                  has to be appended                                  *
*                  >>   real*8    QCHAR1(*) <<                         *
*     NBETA     :  List of residue pointers                  (OUTPUT)  *
*                  >> integer NBETA(NTD,2) <<                          *
*                  NBETA(N,1)=Number of the residue                    *
*                     to which the atom N belongs                      *
*                  NBETA(N,2)=Type of residue                          *
*                     to which the atom N belongs                      *
*     NTD       :  Physical row dimension of NBETA           (OUTPUT)  *
*                  >>   integer NTD <<                                 *
*     NSHIFT    :  Number of atoms currently appended         (INPUT)  *
*                  >> integer NSHIFT <<                                *
*                                                                      *
*---- Last update 25/10/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Universite' de Paris-Sud        *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER natom,nres,atm,ntd,nshift
      INTEGER nbeta(ntd,2)
      CHARACTER*7 label(*),type(*)
      CHARACTER*7 alpha(*),beta(*)
      REAL*8 charg1(*),qchar1(*),mass1(*),qmass1(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n

*==================== EXECUTABLE STATEMENTS ============================


      DO n=1,natom
          alpha(n+nshift)=label(n)
          beta(n+nshift)=type(n)
          qchar1(n+nshift)=charg1(n)
          qmass1(n+nshift)=mass1(n)
          nbeta(n+nshift,1)=atm
          nbeta(n+nshift,2)=nres
      END DO

*-----------------------------------------------------------------------


      RETURN
      END
