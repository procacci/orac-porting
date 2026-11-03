      SUBROUTINE appct(cnncta,cnnctb,na,nb,m,natom,nshift)

************************************************************************
*                                                                      *
*     Append connection table to an existing connection table          *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     CNNCTA  : Connection table to be appended.              (INPUT)  *
*               >> integer CNNCTA(NA,M) <<                             *
*     CNNCTB  : Connection table to which CNNCTA has to                *
*               appended.                                     (IN/OUT) *
*               >> integer CNNCTB(NB,M) <<                             *
*     NA      : Physical row dimension of CNNCTA              (INPUT)  *
*               >> integer NA <<                                       *
*     NB      : Physical row dimension of CNNCTB              (INPUT)  *
*               >> integer NB <<                                       *
*     M       : Physical column dimension of CNNCTA, CNNCTB.  (INPUT)  *
*               >> integer M <<                                        *
*     NATOM   : Number of atoms in CNNCTA.                    (INPUT)  *
*               >> integer NATOM <<                                    *
*     NSHIFT  : Number of Atoms in CNNCTA on input.           (INPUT)  *
*               >> integer NSHIFT <<                                   *
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

      INTEGER      na,nb,m,nshift,natom,
     .             cnncta(na,m),cnnctb(nb,m)
      CHARACTER*80 errmsg

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER j,iatom,coorda,coordb

*==== EXECUTABLE STATEMENTS: ==========================================*

*---- APPEND CONNECTION TABLE TO EXISTING C. T.: ----------------------*

      DO 10 iatom=1,natom
          coorda=cnncta(iatom,1)
          cnnctb(iatom+nshift,1)=coorda
          DO 20 j = 2,1+coorda
              cnnctb(nshift+iatom,j)=cnncta(iatom,j) + nshift
20        CONTINUE
10    CONTINUE


*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
