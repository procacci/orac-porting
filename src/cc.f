      LOGICAL FUNCTION cc(I,J,CONNCT,N)

************************************************************************
*                                                                      *
*     Check connection between atom i and j.                           *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     I       : Number of first atom.                         (INPUT)  *
*               >> integer I <<                                        *
*     J       : Number of second atom.                        (INPUT)  *
*               >> integer J <<                                        *
*     CONNCT  : Connection table.                             (INPUT)  *
*               >> integer CONNCT(N,M) <<                              *
*               CONNCT(I,1) must contain the coordination number       *
*               COORDI of atom no. i, i=1,..,natoms,                 *
*               CONNCT(I,2), ... , CONNCT(i,1+COORDI) the atom numbers *
*               of the neighbors.                                      *
*     N       : Physical column dimension of CONNCT.          (INPUT)  *
*               >> integer N <<                                        *
*     M       : Physical row dimension of CONNCT.             (INPUT)  *
*               >> integer M <<                                        *
*                                                                      *
*---- LAST UPDATE: 03/03/89 -------------------------------------------*
*                                                                      *
*     EXTERNALS: none                                                  *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER i,j,n
      INTEGER connct(n,*)

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER jsrch

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- CHECK I and J: --------------------------------------------------*

      IF(i.LT.1.OR.i.GT.n) THEN
          print *,' CC: Number of first atom out of range'
          CC=.FALSE.
          RETURN
      END IF
      IF(j .LT. 1 .OR. j .GT. n) THEN
          print *,' CC: Number of last atom out of range'
          CC=.FALSE.
          RETURN
      END IF

*---- CHECK CONNECTION OF I and J: ------------------------------------*

      CC    =.FALSE.
      DO 10 jsrch=2,1+connct(i,1)
          IF(connct(i,jsrch).EQ.j) THEN
              CC=.TRUE.
              RETURN
          END IF
10    CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
