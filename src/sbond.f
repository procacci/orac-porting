      SUBROUTINE sbond(connct,n1,natom,ibondl,nbond,n2,iret)

************************************************************************
*                                                                      *
*     Find all bonds from the connection table. A bond is              *
*     defined as a two-atom cluster {a,b}, in which b is               *
*     connected to a. Since each bond occurs in the equivalent         *
*     form a-b and b-a, only the first is kept.                        *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     CONNCT  : Connection Table.                             (INPUT)  *
*               >> integer CONNCT(N1,M1) <<                            *
*               CONNCT(I,1)=Coord. number of atom I,                   *
*                             I=1,..,NATOM;                            *
*               CONNCT(I,J)=Neighbours of atom I,                      *
*                             J=2,..,1+CONNCT(I,1).                    *
*     N1      : Physical row dimension of CONNCT.             (INPUT)  *
*               >> integer N1 <<                                       *
*     M1      : Physical column dimension of CONNCT.          (INPUT)  *
*               >> integer M1 <<                                       *
*     NATOM   : Number of atoms.                              (INPUT)  *
*               >> integer NATOM <<                                    *
*     IBONDL  : List with all bonds.                          (OUTPUT) *
*               >> integer IBONDL(2,N2) <<                             *
*               IBONDL(1..2,I) contains the numbers of atom 1-2        *
*               in bond no. I, I=1,..,NBOND.                           *
*     NBOND   : Number of bonds.                              (OUTPUT) *
*               >> integer NBOND <<                                    *
*     N2      : Physical column dimension of IBONDL.          (INPUT)  *
*               >> integer N2 <<                                       *
*                                                                      *
*---- LAST UPDATE: 03/06/89 -------------------------------------------*
*                                                                      *
*     EXTERNALS:                                                       *
*                                                                      *
*     Written by Gerald Kneller Dept 48B, IBM Kingston 1989            *
*                                                                      *
*     - integer function CC                                            *
*                                                                      *
************************************************************************

      IMPLICIT none

*==== DECLARATIONS: ===================================================*

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER n1,natom,nbond,n2,iret,
     .        connct(n1,*),ibondl(2,*)

*---- EXTERNAL FUNCTIONS: ---------------------------------------------*

      INTEGER  CC
      EXTERNAL CC

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER ibond,i1,i2,a,b,coorda,coordb,ipath,joinab
      LOGICAL ifa
      CHARACTER*80 errmsg

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- INITIALIZATION: -------------------------------------------------*

      DO 10 ibond=1,n2
          ibondl(1,ibond)=0
          ibondl(2,ibond)=0
10    CONTINUE

      nbond=0

*---- FIND ALL POSSIBLE BONDS: ----------------------------------------*

      DO 20 i1=1,natom
          a     =i1
          coorda=connct(a,1)
          DO 30 i2=1,coorda
              b     =connct(a,1+i2)
              IF(a.EQ.b) THEN
                  WRITE(6,1000) a,nbond+1
                  iret=1
                  RETURN
              END IF
              DO 40 ibond=1,nbond
                  IF(a.EQ.ibondl(2,ibond).AND.
     x               b.EQ.ibondl(1,ibond)) THEN
                      GOTO 100
                  END IF
40            CONTINUE
              nbond=nbond + 1
              IF(nbond .GT. n2) THEN
		  errmsg='In SBOND: physical dimension of'//
     x  ' ibondl are insufficient. ABORT!'
                  WRITE(6,'(''ibond = '', i5)') nbond
                  CALL xerror(errmsg,80,1,2)
	      END IF 
              ibondl(1,nbond)=a
              ibondl(2,nbond)=b
100           CONTINUE
30        CONTINUE
20    CONTINUE
      iret=0

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

1000  FORMAT(//' The atoms ',i5,' in bond no. ',i5
     &     ,'  are  I D E N T I C A L '/
     &     /'      F A T A L   E R R O R  in sr sbond '///)
      RETURN
      END
