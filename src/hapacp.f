      SUBROUTINE hapacp(label,n1,natom,cbond,n2,lbond,n3,nohyd,nshfta,
     x                 nshift,iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine append a list of hydrogen bond acceptors         *
*     to a table.                                                      *
*                                                                      *
*     LABEL     :  List of labels for each atoms              (INPUT)  *
*                  >> character*7 LABEL(N1) <<                         *
*     N1        :  Physical dimension of LABEL.                        *
*     NATOM     :  Number of atoms.                           (INPUT)  *
*                  >> integer NATOM <<                                 *
*     CBOND     :  list of hydrogen bond acceptors            (INPUT)  *
*                  >> character*7 CBOND(2,N1) <<                       *
*     N2        :  Physical column dimension of CBOND         (INPUT)  *
*                  >> integer N1 <<                                    *
*     LBOND     :  Table of hydrogen bond acceptors          (OUTPUT)  *
*                  >> integer LBOND <<                                 *
*     N3        :  Physical column dimension of LBOND         (INPUT)  *
*                  >> integer N2 <<                                    *
*     NOHYD     :  Number of hydrogen bond acceptors          (INPUT)  *
*                  >> integer NOHYD <<                                 *
*     NSHIFTA   :  Number of hydrogen bond acceptor in        (INPUT)  *
*                  input to lbond                                      *
*                  >> integer NSHIFTA <<                               *
*     NSHIFT    :  Number of atoms already appended           (INPUT)  *
*                  to the protein tables                               *
*                  >> integer NSHIFT <<                                *
*     IRET      :  Return code.                              (OUTPUT)  *
*                  >> integer IRET <<                                  *
*                  IRET=0 Succesfull return                            *
*                  IRET=1 No match between list of atoms and           *
*                         list of acceptors                            *
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

      INTEGER i,j,k,n
      LOGICAL err,ok

*================= EXECUTABLE STATEMENTS ===============================
      i=1
      err=.false.
10    IF((.NOT.err).AND.(i.LE.nohyd)) THEN
          j=1
20        IF((.NOT.err).AND.(j.LE. 2)) THEN
              n=1
              ok=.false.
30            IF((.NOT.ok).AND.(n.LE.natom)) THEN
                  IF(cbond(j,i).EQ.label(n)) THEN
                      lbond(j,i+nshfta)=n+nshift
                      ok=.true.
                  END IF
                  n=n+1
                  GOTO 30
              END IF
              err=.NOT.ok
              j=j+1
              GOTO 20
          END IF
          i=i+1
          GOTO 10
      END IF

      IF(err) THEN
          iret = 1
          errmsg='  Fatal error in HAPACP : acceptor name NOT FOUND ** '
	  WRITE(6,'(2a7)') (cbond(1,i), cbond(2,i),i=1,nohyd)
      ELSE
          iret = 0
      END IF

      RETURN
      END
