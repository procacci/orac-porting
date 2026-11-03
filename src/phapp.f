      SUBROUTINE phapp(alpha,ntype,lphyd,beta,lpnbd,type,iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine set up an integer list of types for all          *
*     the hydrogen bonded interactions.                                *
*                                                                      *
*     ALPHA         :  List of hydrogen bond types            (INPUT)  *
*                      >> character*7 ALPHA(2,*) <<                    *
*                      ALPHA(1,N)=type of the acceptor.                *
*                      ALPHA(2,N)=type of the donor.                   *
*     NTYPE         :  Physical row and column dimensions of  (INPUT)  *
*                      ALPHA, TYPE and BETA.                           *
*                      >> integer NTYPE <<                             *
*     LPHYD         :  Number of hydrogen bond types          (INPUT)  *
*                      >> integer LPHYD <<                             *
*     BETA          :  List of non-bonded atom types          (INPUT)  *
*                      >> character*7 <<                               *
*     LPNBD         :  Number of non-bonded interaction types (INPUT)  *
*                      >> integer LPNBD <<                             *
*     TYPE          :  Integer list of hydrogen bond types   (OUTPUT)  *
*                      >> integer TYPE(NTYPE,NTYPE) <<                 *
*     IRET          :  Return code                           (OUTPUT)  *
*                      >> integer IRET <<                              *
*                      IRET=0 Normal termination                       *
*                      IRET=1 Fatal error has occured                  *
*     ERRMSG        :  Error message                         (OUTPUT)  *
*                      >> character*80 ERRMSG <<                       *
*                                                                      *
*                                                                      *
*----- Last updated 03/31/89 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*========================= DECLARATIONS ================================

      IMPLICIT none

*-------------------------- ARGUMENTS ----------------------------------

      INTEGER lphyd,lpnbd,iret,ntype
      INTEGER type(ntype,ntype)
      CHARACTER*7 alpha(2,*),beta(*)
      CHARACTER*80 errmsg

*----------------------LOCAL VARIABLES ---------------------------------

      INTEGER m,n,i,j
      CHARACTER*7 char1,char2
      LOGICAL err,ok

*====================== EXECUTABLE STATEMENTS ==========================

      n=1
      err=.false.
10    IF((.NOT.err).AND.(n.LE.lphyd)) THEN
          i=1
          ok=.false.
          char1=alpha(1,n)
          char2=alpha(2,n)
20        IF((.NOT.ok).AND.(i.LE.lpnbd)) THEN
              j=1
30            IF((.NOT.ok).AND.(j.LE.lpnbd)) THEN
                  IF(beta(i).EQ.char1.AND.beta(j).EQ.char2) THEN
                      type(i,j)=n
                      ok=.true.
                  END IF
                  j=j+1
                  GOTO 30
              END IF
              i=i+1
              GOTO 20
          END IF
          err=.NOT.ok
          n=n+1
          GOTO 10
      END IF

*----------------------  Return codes ----------------------------------

      IF(err) THEN
          iret=1
          errmsg='  Fatal error in PHAPP :  acceptor type NOT FOUND **'
      ELSE
          iret=0
          errmsg='   '
      END IF
      RETURN
      END
