      subroutine SCTAB(label,natom,cbonds,n1,nbond,connct,n2,m2)

************************************************************************
*                                                                      *
*     Set up connection table from list of atoms and bonds.            *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     LABEL   : List of labels for each atom.                 (INPUT)  *
*               >> character*7 LABEL(NATOM) <<                         *
*     NATOM   : Number of atoms.                              (INPUT)  *
*               >> integer NATOM <<                                    *
*     CBONDS  : List of bonds.                                (INPUT)  *
*               >> character*7 CBONDS(2,N1) <<                         *
*     N1      : Physical row dimension of CBONDS.             (INPUT)  *
*               >> integer N1 <<                                       *
*     NBOND   : Number of bonds.                              (INPUT)  *
*               >> integer NBOND <<                                    *
*     CONNCT  : Connection Table.                             (OUTPUT) *
*               >> integer CONNCT(N,M) <<                              *
*               CONNCT(I,1) = Coord. number of atom I,                 *
*                             I = 1,..,NATOM;                          *
*               CONNCT(I,J) = Neighbours of atom I,                    *
*                             J = 2,..,1+CONNCT(I,1).                  *
*     N2      : Physical row dimension of CONNCT.             (INPUT)  *
*               >> integer N2 <<                                       *
*     M2      : Physical column dimension of CONNCT.          (INPUT)  *
*               >> integer M2 <<                                       *
*                                                                      *
*---- LAST UPDATE: 03/27/89 -------------------------------------------*
*                                                                      *
*     Written by Gerald Kneller Dept 48B, IBM Kingston 1989            *
*                                                                      *
*     EXTERNALS: none                                                  *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      IMPLICIT none

*---- ARGUMENTS: ------------------------------------------------------*

      integer     natom,nbond,n1,n2,m2,
     .            connct(n2,m2)
      character*7 label(natom),cbonds(2,n1)
      CHARACTER*80 errmsg

*---- LOCAL VARIABLES: ------------------------------------------------*

      integer i,j,ibond,iatom,coordi,coordj

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- INITIALISE CONNECTION TABLE: ------------------------------------*

      DO 10 i=1,n2
          DO 20 j=1,m2
              connct(i,j) = 0
20        CONTINUE
10    CONTINUE

*---- SET UP CONNECTION TABLE: ----------------------------------------*

      DO 30 ibond=1,nbond
          i=0
          j=0
          DO 40 iatom = 1,natom
            if (cbonds(1,ibond) .EQ. label(iatom)) i = iatom
            if (cbonds(2,ibond) .EQ. label(iatom)) j = iatom
40        CONTINUE
          IF(i.ne.0.AND.j.ne.0) THEN
             coordi             = connct(i,1) + 1
             IF(coordi .GT. m2-1) THEN
                 errmsg=' Error in SCTAB: connection tables exceed '//
     x                  'physical dimensions. Abort.'
                 CALL xerror(errmsg,80,1,2)
             END IF
             connct(i,1)        = coordi
             connct(i,1+coordi) = j
             coordj             = connct(j,1) + 1
             IF(coordj .GT. m2-1) THEN
                 errmsg=' Error in SCTAB: connection tables exceed '//
     x                  'physical dimensions. Abort.'
                 CALL xerror(errmsg,80,1,2)
             END IF
             connct(j,1)        = coordj
             connct(j,1+coordj) = i
         ELSE
             errmsg=' Error in SCTAB: Bond label not in atom list.'//
     x ' Abort. '
             CALL xerror(errmsg,80,1,2)
         END IF
30    CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
