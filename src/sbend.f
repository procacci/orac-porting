      SUBROUTINE sbend(connct,n1,natom,ibendl,nbend,n2)

************************************************************************
*                                                                      *
*     Find all bendings from the connection table. A bending is        *
*     defined as a three-atom cluster {a,b,c}, in which b is           *
*     connected to a and c connected to b. Since each bending          *
*     occurs in the equivalent forms a-b-c and c-b-a, only the the     *
*     first is kept. The bending angle is the angle between the        *
*     connections a-b and b-c.                                         *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     CONNCT  : Connection Table.                             (INPUT)  *
*               >> integer CONNCT(N1,M1) <<                            *
*               CONNCT(I,1)=Coord. number of atom I,                 *
*                             I=1,..,NATOM;                          *
*               CONNCT(I,J)=Neighbours of atom I,                    *
*                             J=2,..,1+CONNCT(I,1).                  *
*     N1      : Physical row dimension of CONNCT.             (INPUT)  *
*               >> integer N1 <<                                       *
*     M1      : Physical column dimension of CONNCT.          (INPUT)  *
*               >> integer M1 <<                                       *
*     NATOM   : Number of atoms.                              (INPUT)  *
*               >> integer NATOM <<                                    *
*     IBENDL  : List with all bendings.                       (OUTPUT) *
*               >> integer IBENDL(3,N2) <<                             *
*               IBENDL(1..3,I) contains the numbers of atom 1-2-3      *
*               in bending no. I, I=1,..,NBEND.                      *
*     NBEND   : Number of bends.                              (OUTPUT) *
*               >> integer NBEND <<                                    *
*     N2      : Physical row dimension of IBENDL.             (INPUT)  *
*               >> integer N2 <<                                       *
*                                                                      *
*---- LAST UPDATE: 03/06/89 -------------------------------------------*
*                                                                      *
*     Written by Gerald Kneller Dept 48B, IBM Kingston 1989            *
*                                                                      *
*     EXTERNALS:                                                       *
*                                                                      *
*     - integer function CC                                            *
*                                                                      *
************************************************************************

      IMPLICIT none

*==== DECLARATIONS: ===================================================*

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER n1,n2,natom,nbend,
     .        connct(n1,*),ibendl(3,*)

*---- EXTERNAL FUNCTIONS: ---------------------------------------------*

      INTEGER  CC
      EXTERNAL CC

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER ibend,i1,i2,i3,i4,a,b,c,d,coorda,coordb,coordc
      LOGICAL alldif
      CHARACTER*80 errmsg

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- INITIALIZATION: -------------------------------------------------*

      DO 10 ibend=1,n2
          ibendl(1,ibend)=0
          ibendl(2,ibend)=0
          ibendl(3,ibend)=0
10    CONTINUE

      nbend=0

*---- FIND ALL POSSIBLE BENDINGS: -------------------------------------*

      DO 20 i1=1,natom
          a     =i1
          coorda=connct(a,1)
          DO 30 i2=1,coorda
              b     =connct(a,1+i2)
              coordb=connct(b,1)
              DO 40 i3=1,coordb
                  c     =connct(b,1+i3)
                  IF(a .ne. b .AND. a .ne. c  .AND. b .ne. c) THEN
                      coordc=connct(c,1)
                      DO 60 i4=1,coordc
                          d     =connct(c,1+i4)
                          IF(d .EQ. a) GOTO 100
60                    CONTINUE
                      alldif=.true.
                  ELSE
                      alldif=.false.
                  END IF
                  IF(alldif) THEN
                      DO 50 ibend=1,nbend
                          if(a .EQ. ibendl(3,ibend) .AND.
     .                       b .EQ. ibendl(2,ibend) .AND.
     .                       c .EQ. ibendl(1,ibend)      ) THEN
                              GOTO 100
                          END IF
50                    CONTINUE
                      nbend=nbend + 1
		      IF(nbend .GT. n2) THEN
			  errmsg='In SBEND: physical dimension of'//
     x  ' ibendl are insufficient. ABORT!'
                          CALL xerror(errmsg,80,1,2)
	              END IF 
                      ibendl(1,nbend)=a
                      ibendl(2,nbend)=b
                      ibendl(3,nbend)=c
                  END IF
100               CONTINUE
40            CONTINUE
30        CONTINUE
20    CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
