      SUBROUTINE pbend(ibendl,nbend,n1,type,tbendp,pbendp,nbendp,n2,
     x                 pbendx,m2,bug,nbug)

************************************************************************
*                                                                      *
*     Check bendings with potential and generate parameter list for    *
*     bendings.                                                        *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     IBENDL  : List with all bendings found from sr. SBEND.  (INPUT)  *
*               >> integer IBENDL(3,N1) <<                             *
*               IBENDL(1..3,I) contain the numbers of atom 1-2-3       *
*               in bending no. I, I = 1,..,NBEND.                      *
*     NBEND   : Number of bends.                              (INPUT)  *
*               >> integer NBEND <<                                    *
*     N1      : Physical column dim. IBENDL/PBENDX.           (INPUT)  *
*               >> integer N1 <<                                       *
*     TYPE    : Types for each atom in the bending list.      (INPUT)  *
*               >> CHARACTER*7 TYPE(*) <<                              *
*     TBENDP  : List of bendings in the potential model.      (INPUT)  *
*               >> character*7 TBENDP(3,N2) <<                         *
*               TBENDP(1..3,I) contains the types of atom 1-2-3        *
*               in potential no. I, I = 1,..,NBENDP.                   *
*     PBENDP  : List with the potential parameters            (INPUT)  *
*               >> real PBENDP(N2,M2) <<                               *
*     NBENDP  : Number of potential parameters                (INPUT)  *
*               >> integer NBENDP <<                                   *
*     N2      : Physical row dimension of PBENDP and column            *
*               dimension of TBENDP                                    *
*               >> integer N2 <<                                       *
*     PBENDX  : List with the potential parameters            (OUTPUT) *
*               corresponding to the bendings in IBENDL.               *
*               >> real PBENDK(M2,N1) <<                               *
*     M2      : Physical column dimension of PBENDP/PBENDX    (INPUT)  *
*               >> integer M2 <<                                       *
*     IRET    : Return code.                                  (OUTPUT) *
*               >> integer IRET <<                                     *
*               IRET = 0,..,3. For each value of IRET a                *
*               corresponding error message is assigned to             *
*               ERRMSG (see source code).                              *
*     ERRMSG  : Error message.                                (OUTPUT) *
*               >> CHARACTER*80 ERRMSG <<                              *
*                                                                      *
*---- LAST UPDATE: 03/30/89 -------------------------------------------*
*                                                                      *
*                                                                      *
*     Written by Gerald Kneller Dept 48B, IBM Kingston 1989            *
*                                                                      *
*     EXTERNALS:                                                       *
*                                                                      *
*     - NONE                                                           *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      IMPLICIT CHARACTER*80(a-z)

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER      nbend,nbendp,n1,n2,m2,iret,ibendl(3,n1),bug(*),nbug
      REAL*8       pbendx(n1,*),pbendp(n2,*)
      CHARACTER*80 errmsg
      CHARACTER*7  type(*),tbendp(3,n2)

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER      ibend,ibendp,nbendk,m,j1,j2

*==== EXECUTABLE STATEMENTS: ==========================================*



*---- INITIALIZATION:  ------------------------------------------------*

      DO 10 ibend = 1,n1
          DO 20 j2 = 1,m2
              pbendx(ibend,j2) = 0.0D0
20        CONTINUE
10    CONTINUE


      nbug=0
      nbendk = 0
      DO 30 ibend = 1,nbend
          DO 40 ibendp = 1,nbendp
              IF(type(ibendl(1,ibend)) .EQ. tbendp(1,ibendp) .AND.
     x          type(ibendl(2,ibend)) .EQ. tbendp(2,ibendp) .AND.
     x          type(ibendl(3,ibend)) .EQ. tbendp(3,ibendp)     .OR.
     x          type(ibendl(1,ibend)) .EQ. tbendp(3,ibendp) .AND.
     x          type(ibendl(2,ibend)) .EQ. tbendp(2,ibendp) .AND.
     x          type(ibendl(3,ibend)) .EQ. tbendp(1,ibendp)) THEN
                  nbendk = nbendk + 1
                  DO 50 m=1,m2
                      pbendx(nbendk,m) = pbendp(ibendp,m)
50                CONTINUE
                  GOTO 100
              END IF
40        CONTINUE
          nbug=nbug+1
          bug(nbug)=ibend 
c---      one thousend Bugs: we had enough 
          if(nbug.gt.1000) RETURN 
100       CONTINUE
30    CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN

      END 

