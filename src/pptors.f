      SUBROUTINE pptors(itorsl,ntors,o1,type,ttorsp,ptorsp,ntorsp,o2,
     x                 ptorsx,p2,bug,nbug)

************************************************************************
*                                                                      *
*     Check p-torsions with potential and generate parameter list for  *
*     torsions.                                                        *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     ITORSL  : List with all p-torsions found from sr. STORS.         *
*               >> integer ITORSL(4,N1) <<                             *
*               ITORSL(1..4,I) contain the numbers of atom 1-2-3-4     *
*               in torsion no. I, I = 1,..,NTORS.    (INPUT/OUTPUT)    *
*                                                                      *
*     NTORS   : Number of torsions.                  (INPUT/OUTPUT)    *
*               >> integer NTORS <<                                    *
*     N1      : Physical column dim. ITORSL/PTORSX.           (INPUT)  *
*               >> integer N1 <<                                       *
*     TYPE    : Types for each atom in the torsions list.     (INPUT)  *
*               >> CHARACTER*7 TYPE(*) <<                              *
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
*---- LAST UPDATE: 06/28/89 -------------------------------------------*
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
      REAL*8       ptorsx(o1,p2),ptorsp(o2,p2),eps
      CHARACTER*80 errmsg
      CHARACTER*7  match1(4),match2(4),type(*),ttorsp(4,o2)
      CHARACTER*1 wild

*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER na,work
      PARAMETER (na=m3)
      COMMON /rag1/ work(4,na)
      INTEGER      itors,itorsp,ntorsk,m,j,j1,j2
      LOGICAL    err
      DATA wild/'x'/
      DATA eps/1.0D-5/

*==== EXECUTABLE STATEMENTS: ==========================================*


*---- INITIALIZATION:  ------------------------------------------------*

      iret   = 0
      DO 10 itors = 1,o1
          DO 20 j2 = 1,p2
              ptorsx(itors,j2) = 0.0D0
20        CONTINUE
10    CONTINUE

*=======================================================================

*=======================================================================
*------- Look for all the matchin torsions -----------------------------
*=======================================================================

      ntorsk = 0
      nbug=0
      DO 30 itors = 1,ntors

*=======================================================================
*------- Keep matching torsions without wildcards ----------------------
*------- Do not retain those torsions with a null potential ------------
*=======================================================================

          err=.TRUE.
          DO 40 itorsp = 1,ntorsp
              IF(type(itorsl(1,itors)) .EQ. ttorsp(1,itorsp) .AND.
     x            type(itorsl(2,itors)) .EQ. ttorsp(2,itorsp) .AND.
     x            type(itorsl(3,itors)) .EQ. ttorsp(3,itorsp) .AND.
     x            type(itorsl(4,itors)) .EQ. ttorsp(4,itorsp) .OR.
     x            type(itorsl(1,itors)) .EQ. ttorsp(4,itorsp) .AND.
     x            type(itorsl(2,itors)) .EQ. ttorsp(3,itorsp) .AND.
     x            type(itorsl(3,itors)) .EQ. ttorsp(2,itorsp) .AND.
     x            type(itorsl(4,itors)) .EQ. ttorsp(1,itorsp)) THEN
                  ntorsk = ntorsk + 1
                  IF(ntorsk .GT. na) THEN
                      errmsg='In pptors: Physical dimensions of the'
     x                //' work array are insufficient. '
                      iret=1
                      RETURN
                  END IF
                  DO 50 m=1,p2
                      ptorsx(ntorsk,m) = ptorsp(itorsp,m)
50                CONTINUE
                  DO 60 m=1,4
                      work(m,ntorsk) = itorsl(m,itors)
60                CONTINUE
                  err=.FALSE.
              END IF
40        CONTINUE
          IF(.NOT.err) GOTO 100

*=======================================================================
*----------- Keep matching torsions with wildcard ----------------------
*------- Do not retain those torsions with a null potential ------------
*=======================================================================

          DO 70 itorsp=1,ntorsp
              DO 80 j=1,4
                  IF(ttorsp(j,itorsp).EQ.wild) THEN
                      match1(j)=type(itorsl(j,itors))
                  ELSE
                      match1(j)=ttorsp(j,itorsp)
                  END IF
80            CONTINUE
              DO 90 j=1,4
                  IF(ttorsp(5-j,itorsp).EQ.wild) THEN
                      match2(j)=type(itorsl(j,itors))
                  ELSE
                      match2(j)=ttorsp(5-j,itorsp)
                  END IF
90            CONTINUE
              IF(type(itorsl(1,itors)) .EQ. match1(1) .AND.
     x          type(itorsl(2,itors)) .EQ. match1(2) .AND.
     x          type(itorsl(3,itors)) .EQ. match1(3) .AND.
     x          type(itorsl(4,itors)) .EQ. match1(4) .OR.
     x          type(itorsl(1,itors)) .EQ. match2(1) .AND.
     x          type(itorsl(2,itors)) .EQ. match2(2) .AND.
     x          type(itorsl(3,itors)) .EQ. match2(3) .AND.
     x          type(itorsl(4,itors)) .EQ. match2(4)) THEN
!                 IF(DABS(ptorsp(itorsp,1)).GT.eps) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na*3) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF
                      DO 110 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
110                   CONTINUE
                      DO 120 m=1,4
                          work(m,ntorsk) = itorsl(m,itors)
120                   CONTINUE
!                 END IF
                  GOTO 100
              END IF
70        CONTINUE
          nbug=nbug+1
          bug(nbug)=itors
c---      one thousend Bugs: we had enough 
          if(nbug.gt.1000) RETURN 
100       CONTINUE
30    CONTINUE

      IF(nbug.gE.1) RETURN

*=======================================================================
*---------- Change the list of p-torsions ------------------------------
*=======================================================================

      ntors=ntorsk
      DO 130 j=1,ntors
          DO 140 j1=1,4
              itorsl(j1,j)=work(j1,j)
140       CONTINUE
130   CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN

      END
