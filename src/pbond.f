      SUBROUTINE pbond(ibondl,nbond,n1,type,tbondp,pbondp,nbondp,n2,
     x                 pbondx,m2,ibond_bug,nbondk_bug)

************************************************************************
*                                                                      *
*     Check bondings with potential and generate parameter list for    *
*     bondings.                                                        *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     IBONDL  : List with all bondings found from sr. SBOND.  (INPUT)  *
*               >> integer IBONDL(3,N1) <<                             *
*               IBONDL(1..3,I) contain the numbers of atom 1-2-3       *
*               in bonding no. I, I = 1,..,NBOND.                      *
*     NBOND   : Number of bonds.                              (INPUT)  *
*               >> integer NBOND <<                                    *
*     N1      : Physical column dim. IBONDL/PBONDX.           (INPUT)  *
*               >> integer N1 <<                                       *
*     TYPE    : Types for each atom in the bonding list.      (INPUT)  *
*               >> CHARACTER*7 TYPE(*) <<                              *
*     TBONDP  : List of bondings in the potential model.      (INPUT)  *
*               >> character*7 TBONDP(3,N2) <<                         *
*               TBONDP(1..3,I) contains the types of atom 1-2-3        *
*               in potential no. I, I = 1,..,NBONDP.                   *
*     PBONDP  : List with the potential parameters            (INPUT)  *
*               >> real PBONDP(N2,M2) <<                               *
*     NBONDP  : Number of potential parameters                (INPUT)  *
*               >> integer NBONDP <<                                   *
*     N2      : Physical row dimension of PBONDP and column            *
*               dimension of TBONDP                                    *
*               >> integer N2 <<                                       *
*     PBONDX  : List with the potential parameters            (OUTPUT) *
*               corresponding to the bondings in IBONDL.               *
*               >> real PBONDK(M2,N1) <<                               *
*     M2      : Physical column dimension of PBONDP/PBONDX    (INPUT)  *
*               >> integer M2 <<                                       *
*                                                                      *
*---- LAST UPDATE: 05/17/94 -------------------------------------------*
*                                                                      *
*                                                                      *
*     Written by Massimo Marchi CEA, CE Saclay FRANCE 1994             *
*                                                                      *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      IMPLICIT none

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER      nbond,nbondp,n1,n2,m2,ibondl(2,*),ibond_bug(*)
     &     ,nbondk_bug
      REAL*8       pbondx(n1,*),pbondp(n2,*)
      CHARACTER*7  type(*),tbondp(2,*)

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER      ibond,ibondp,nbondk,m,j2

*==== EXECUTABLE STATEMENTS: ==========================================*



*---- INITIALIZATION:  ------------------------------------------------*

      DO 10 ibond = 1,n1
          DO 20 j2 = 1,m2
              pbondx(ibond,j2) = 0.0D0
20        CONTINUE
10    CONTINUE


      nbondk = 0
      nbondk_bug = 0
      DO ibond = 1,nbond
          DO ibondp = 1,nbondp
              IF(type(ibondl(1,ibond)) .EQ. tbondp(1,ibondp) .AND.
     x          type(ibondl(2,ibond)) .EQ. tbondp(2,ibondp) .OR.
     x          type(ibondl(1,ibond)) .EQ. tbondp(2,ibondp) .AND.
     x          type(ibondl(2,ibond)) .EQ. tbondp(1,ibondp)) THEN
                  nbondk = nbondk + 1
                  DO m=1,m2
                      pbondx(nbondk,m) = pbondp(ibondp,m)
                  END DO
                  GOTO 100
              END IF
          END DO
          nbondk_bug=nbondk_bug+1
          ibond_bug(nbondk_bug)=ibond 
c---      one thousand Bugs: we had enough 
          if(nbondk_bug.gt.1000) RETURN 
100       CONTINUE
      END DO

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
