      SUBROUTINE hapdon(label,n1,natom,cbond,n2,lbond,n3,nohyd,nshfta,
     x                 nshift,iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine appends a list of hydrogen bond donors           *
*     to a table.                                                      *
*                                                                      *
*     LABEL     :  List of labels for each atoms              (INPUT)  *
*                  >> character*7 LABEL(N1) <<                         *
*     N1        :  Physical dimension of LABEL.                        *
*                  >> integer N1 <<                                    *
*     NATOM     :  Number of atoms.                           (INPUT)  *
*                  >> integer NATOM <<                                 *
*     CBOND     :  list of hydrogen bond donors               (INPUT)  *
*                  >> character*7 CBOND(2,N2) <<                       *
*     N2        :  Physical column dimension of CBOND         (INPUT)  *
*                  >> integer N2 <<                                    *
*     LBOND     :  Table of hydrogen bond donors             (OUTPUT)  *
*                  >> integer LBOND(2,N3) <<                           *
*     N3        :  Physical column dimension of LBOND         (INPUT)  *
*                  >> integer N3 <<                                    *
*     NOHYD     :  Number of hydrogen bond donors             (INPUT)  *
*                  >> integer NOHYD <<                                 *
*     NSHIFTA   :  Number of hydrogen bond donor in           (INPUT)  *
*                  input to lbond                                      *
*                  >> integer NSHIFTA <<                               *
*     NSHIFT    :  Number of atoms already appended           (INPUT)  *
*                  to the protein tables                               *
*                  >> integer NSHIFT <<                                *
*     IRET      :  Return code.                              (OUTPUT)  *
*                  >> integer IRET <<                                  *
*                  IRET=0 Succesfull return                            *
*                  IRET=1 No match between list of atoms and           *
*                         list of donors                               *
*     ERRMSG    :  Error message                             (OUTPUT)  *
*                 >> character*80 ERRMSG <<                            *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*=================== DECLARATIONS ======================================

      IMPLICIT CHARACTER*80(a-z)

*--------------- Arguments ---------------------------------------------

      INTEGER natom,n1,n2,n3,nohyd,nshfta,nshift,iret,lbond(2,n3)
      CHARACTER*7 label(n1),cbond(2,n2)
      CHARACTER*80 errmsg

*--------------- Local arguments ---------------------------------------

      INTEGER i,j,n

*================= EXECUTABLE STATEMENTS ===============================

      DO 10 i=1,nohyd
          DO 20 j=1,2

              DO 30 n=1,natom
                  IF(cbond(j,i).EQ.label(n)) THEN
                      lbond(3-j,i+nshfta)=n+nshift
                      GOTO 100
                  END IF
30            CONTINUE
              GOTO 200
100           CONTINUE
20        CONTINUE
10    CONTINUE

*---------------- Succesfully return to the calling routine ------------

      iret=0
      RETURN

*--------------- Return code 1 -----------------------------------------

200   CONTINUE
      iret=1
      errmsg='  Fatal error in HAPDON : donor    name NOT FOUND ** '
      WRITE(6,'(2a7)') (cbond(1,i), cbond(2,i),i=1,nohyd)
      RETURN
      END
