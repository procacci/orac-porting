      SUBROUTINE pitors(itorsl,ntors,o1,type,ttorsp,ptorsp,ntorsp,o2,
     x                 ptorsx,p2,bug,nbug)

************************************************************************
*                                                                      *
*     Check improper torsions with potential and generate parameter    *
*     list for improper torsions.                                      *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     ITORSL  : List with all p-torsions found from sr. SITORS.        *
*               >> integer ITORSL(4,N1) <<                             *
*               ITORSL(1..4,I) contain the numbers of atom 1-2-3-4     *
*               in torsion no. I, I = 1,..,NTORS.             (INPUT)  *
*                                                                      *
*     NTORS   : Number of torsions.                           (INPUT)  *
*               >> integer NTORS <<                                    *
*     N1      : Physical column dim. ITORSL/PTORSX.           (INPUT)  *
*               >> integer N1 <<                                       *
*     TYPE    : Types for each atom in the torsions list.     (INPUT)  *
*               >> character*7 TYPE(*) <<                              *
*     TTORSP  : List of torsions in the potential model.      (INPUT)  *
*               >> character*7 TTORSP(4,N2) <<                         *
*               TTORSP(1..4,I) contains the types of atom 1-2-3-4      *
*               in potential no. I, I = 1,..,NTORSP.                   *
*     PTORSP  : List with the potential parameters            (INPUT)  *
*               >> real PTORSP(N2,M2) <<                               *
*     NTORSP  : Number of potential parameters                (INPUT)  *
*               >> integer NTORSP <<                                   *
*     N2      : Physical row dimension of PTORSP and column            *
*               dimension of TTORSP                                    *
*               >> integer N2 <<                                       *
*     PTORSX  : List with the potential parameters            (OUTPUT) *
*               corresponding to the torsions in ITORSL.               *
*               >> real PTORSK(M2,N1) <<                               *
*     M2      : Physical column dimension of PTORSP/PTORSX    (INPUT)  *
*               >> integer M2 <<                                       *
*     WORK    : Work array                                    (INPUT)  *
*               >> integer WORK(4,N1) <<                               *
*               In the calling program must have the same              *
*               physical column dimension as ITORSL                    *
*     IRET    : Return code.                                  (OUTPUT) *
*               >> integer IRET <<                                     *
*               IRET = 0,..,1. For each value of IRET a                *
*               corresponding error message is assigned to             *
*               ERRMSG.                                                *
*     ERRMSG  : Error message.                                (OUTPUT) *
*               >> CHARACTER*80 ERRMSG <<                              *
*                                                                      *
*---- LAST UPDATE: 05/24/89 -------------------------------------------*
*                                                                      *
*     Modified version of a subroutine written by                      *
*     Gerald Kneller Dept 48B, IBM Kingston 1989                       *
*                                                                      *
*     EXTERNALS:                                                       *
*                                                                      *
*     - NONE                                                           *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      use parst

      IMPLICIT none

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER      ntors,ntorsp,o1,o2,p2,iret,itorsl(4,o1),bug(*),nbug
      REAL*8       ptorsx(o1,p2),ptorsp(o2,p2)
      CHARACTER*80 errmsg
      CHARACTER*7  match1(4),type(*),ttorsp(4,o2)
      CHARACTER*1 wild

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER na,work
      PARAMETER (na=m4)
      COMMON /rag1/ work(4,na)
      INTEGER      itors,itorsp,ntorsk,m,j,j1,j2
      INTEGER      ip,iperm(6,4)
      LOGICAL err
      DATA iperm     /1,2,4,4,1,2,
     X                2,1,1,2,4,4,
     X                3,3,3,3,3,3,
     X                4,4,2,1,2,1/
      DATA wild/'x'/
      character*4 tt1,tt2,tt3,tt4

*==== EXECUTABLE STATEMENTS: ==========================================*


*---- INITIALIZATION:  ------------------------------------------------*nn

      DO 10 itors = 1,o1
          DO 20 j2 = 1,p2
              ptorsx(itors,j2) = 0.0D0
20        CONTINUE
10    CONTINUE

*=======================================================================

*------- Look for all the matching torsions ----------------------------

      nbug = 0
      ntorsk = 0
      DO 30 itors = 1,ntors

*------- Keep matching improper torsions -------------------------------

          DO 90 itorsp = 1,ntorsp
              DO 110 ip=1,6
                  DO 120 j=1,4
                      IF(ttorsp(j,itorsp).EQ.wild) THEN
                          match1(j)=type(itorsl(iperm(ip,j),itors))
                      ELSE
                          match1(j)=ttorsp(j,itorsp)
                      END IF
120               CONTINUE
                  IF( type(itorsl(iperm(ip,1),itors)).EQ.match1(1).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(4)) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na*3) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF
                      DO 130 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
130                   CONTINUE
                      DO 140 m=1,4
                          work(m,ntorsk) = itorsl(m,itors)
140                   CONTINUE
                      GOTO 100
                  END IF
110           CONTINUE
90        CONTINUE
          nbug=nbug+1
          bug(nbug)=itors
c---      one thousend Bugs: we had enough 
          if(nbug.gt.1000) RETURN 
100       CONTINUE
30    CONTINUE

      if(nbug.ge.1) RETURN

*---------- Change the list of p-torsions ------------------------------

      ntors=ntorsk
      DO 250 j=1,ntors
          DO 260 j1=1,4
              itorsl(j1,j)=work(j1,j)
260       CONTINUE
250   CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
