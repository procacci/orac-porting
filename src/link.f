      SUBROUTINE linka(label,nato,clinks,n1,nlink,connct,n2,
     x                resid)

************************************************************************
*                                                                      *
*     Link submolecules to supermolecule.                              *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     LABEL   : List of labels for each atom in the                    *
*               supermolecule.                                (INPUT)  *
*               >> character*7  LABEL(NATOM) <<                        *
*     NATO    : Number of Atomes                              (INPUT)  *
*     CLINKS  : List of links.                                (INPUT)  *
*               >> character*16 CLINKS(N1,2) <<                        *
*     N1      : Physical row dimension of CLINKS.             (INPUT)  *
*               >> integer N1 <<                                       *
*     NLINK   : Number of links, not includin extra links.    (INPUT)  *
*               >> integer NLINK <<                                    *
*     CONNCT  : Connection table for the supermolecule.       (IN/OUT) *
*               >> integer CONNCT(N,M) <<                              *
*               CONNCT(I,1)=Coord. number of atom I,                   *
*                             I=1,..,NATOM;                            *
*               CONNCT(I,J)=Neighbours of atom I,                      *
*                             J=2,..,1+CONNCT(I,1).                    *
*     N2      : Physical row dimension of CONNCT.             (INPUT)  *
*               >> integer N2 <<                                       *
*     RESID   : Residue identifier.                                    *
*               >> integer RESID(NATOM) <<                        (I)  *
*                                                                      *
*---- LAST UPDATE: 05/18/89 -------- M. Marchi ------------------------*
*                                                                      *
*     Adapted  version of a subroutine written by                      *
*     Gerald Kneller Dept 48B, IBM Kingston 1989                       *
*                                                                      *
*     EXTERNALS: none                                                  *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      IMPLICIT none

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER      nato,nlink,n1,n2,
     .             connct(n2,*),resid(*)
      CHARACTER*7 label(*),clinks(n1,*)

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER i,j,ilink,iatom,coordi,coordj,resp,resn,iares

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- SET LINKS IN CONNECTION TABLE: ----------------------------------*

      DO 10 ilink=1,nlink
          i=0
          j=0
          resp=ilink
          resn=ilink+1
          DO 20 iatom=1,nato
              iares=resid(iatom)
              IF(clinks(ilink,1).EQ.label(iatom).AND.
     x           iares.EQ.resp) i=iatom
              IF(clinks(ilink,2).EQ.label(iatom).AND.
     x           iares.EQ.resn) j=iatom
20        CONTINUE
          IF(i .ne. 0 .AND. j .ne. 0) THEN
              coordi            =connct(i,1) + 1
              connct(i,1)       =coordi
              connct(i,1+coordi)=j
              coordj            =connct(j,1) + 1
              connct(j,1)       =coordj
              connct(j,1+coordj)=i
          END IF
10    CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
