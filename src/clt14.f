      SUBROUTINE clt14(connct,o1,p1,natom,lbndg,lbend,int14,int14p,o2,
     x                 itorsl,ntors,nprot,protl,dihed)

************************************************************************
*                                                                      *
*     CLT14 compute a list of 1-4 interactions.                        *
*     Find all proper torsions from the connection table. A proper     *
*     torsion is defined as a four atom cluster {a,b,c,d}, in which    *
*     b is connected to a, c connected to b, and d connected to c.     *
*     The torsion angle is the angle between the planes a-b-c and      *
*     b-c-d. Since each torsion occurs in the two equivalent forms     *
*     a-b-c-d and d-c-b-a, only the first is kept.                     *
*                                                                      *
*                                                                      *
*     CONNCT  : Connection Table.                             (INPUT)  *
*               >> INTEGER CONNCT(N1,M1) <<                            *
*               CONNCT(I,1)=Coord. number of atom I,                   *
*                             I=1,..,NATOM;                            *
*               CONNCT(I,J)=Neighbours of atom I,                      *
*                             J=2,..,1+CONNCT(I,1).                    *
*     N1      : Physical row dimension of CONNCT.             (INPUT)  *
*               >> INTEGER N1 <<                                       *
*     M1      : Physical column dimension of CONNCT.          (INPUT)  *
*               >> INTEGER M1 <<                                       *
*     NATOM   : Number of atoms.                              (INPUT)  *
*               >> INTEGER NATOM <<                                    *
*     LBNDG   :  List of all bendings.                            (I)  *
*                >> integer LBNDG(3,*) <<                              *
*     LBEND   :  Number of bendings.                              (I)  *
*     INT14   :  List of 1-4 interactions.                        (O)  *
*             :  >> integer INT14(2,*) <<                              *
*     INT14P  :  Number of 1-4 interactions.                      (O)  *
*     N2      : Physical row dimension of ITORSL.             (INPUT)  *
*               >> INTEGER N2 <<                                       *
*                                                                      *
*----- Last update 02/21/91 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley                            *
*                                                                      *
*     EXTERNALS none                                                   *
*                                                                      *
*                                                                      *
************************************************************************

      IMPLICIT none

*==== DECLARATIONS: ===================================================*

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER o1,p1,natom,ntors,o2,protl(*),nprot,
     .        connct(o1,p1)
      INTEGER lbend,int14p
      INTEGER lbndg(3,*),int14(2,*),itorsl(4,*)
      LOGICAL dihed,la4,ld1,la1,ld4,lb2,lc3,lb3,lc2,lad14,lbc23

*---- WORK COMMON BLOCK -----------------------------------------------*


*---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER itors,i1,i2,i3,i4,a,b,c,d,coorda,coordb,coordc,
     .        ipath,e,f,g,coorde,coordf
      INTEGER i,j,n,count,offset,cprot,ii
      LOGICAL alldif
      LOGICAL ok
      CHARACTER*80 errmsg


*==== EXECUTABLE STATEMENTS: ==========================================*

*---- INITIALIZATION: -------------------------------------------------*

      DO 10 itors=1,o2
          itorsl(1,itors)=0
          itorsl(2,itors)=0
          itorsl(3,itors)=0
          itorsl(4,itors)=0
10    CONTINUE

      ntors=0

*---- FIND ALL POSSIBLE PROPER TORSIONS: ------------------------------*

      offset=0
      DO i=1,nprot
         cprot=protl(1+offset)
         count=0
         DO 20 ii=1,cprot
            i1=protl(1+offset+ii)
            a     =i1
            coorda= connct(a,1)
            DO 30 i2=1,coorda
               b     =connct(a,1+i2)
               coordb=connct(b,1)
               DO 40 i3=1,coordb
                  c     =connct(b,1+i3)
                  coordc=connct(c,1)
                  IF(a .ne. c) THEN
                     DO 50 i4=1,coordc
                        d    =connct(c,1+i4)
                        IF(d .NE. b .AND. a .NE. d) THEN
                           alldif=.true.
                        ELSE
                           alldif=.false.
                        END IF
                        IF(alldif) THEN
                           DO 60 itors=ntors-count+1,ntors
                              la4 =  a .EQ. itorsl(4,itors)
                              ld1 =  d .EQ. itorsl(1,itors)
                              la1 =  a .EQ. itorsl(1,itors)
                              ld4 =  d .EQ. itorsl(4,itors)
                              lb2 =  b .EQ. itorsl(2,itors)
                              lc3 =  c .EQ. itorsl(3,itors)
                              lb3 =  b .EQ. itorsl(3,itors)
                              lc2 =  c .EQ. itorsl(2,itors)
                              lad14 = (la4.and.ld1).or.(la1.and.ld4)
                              lbc23 = (lb2.and.lc3).or.(lb3.and.lc2)
                              if(lad14.and.lbc23) goto 100
60                         CONTINUE
                           ntors=ntors + 1
                           count=count+1
                           IF(o2 .LT. ntors) THEN
                              errmsg='In CLT14: Physical dimensions'//
     x' of the work array are insufficient. ABORT'
                              CALL xerror(errmsg,80,1,2)
                           END IF
                           itorsl(1,ntors)=a
                           itorsl(2,ntors)=b
                           itorsl(3,ntors)=c
                           itorsl(4,ntors)=d
                        END IF
100                     CONTINUE
50                   CONTINUE
                  END IF
40             CONTINUE
30          CONTINUE
20       CONTINUE
         offset=offset+cprot+1
      END DO

      n=0
      DO i=1,ntors
          i1=itorsl(1,i)
          i2=itorsl(4,i)
          ok=.TRUE.
          DO j=1,lbend
              IF((lbndg(1,j).EQ.i1.AND.lbndg(3,j).EQ.i2).OR.
     x           (lbndg(3,j).EQ.i1.AND.lbndg(1,j).EQ.i2)) THEN
                  ok=.FALSE.
              END IF
              IF (lbndg(1,j).EQ.i1.AND.lbndg(2,j).EQ.i2) ok=.false.
              IF (lbndg(2,j).EQ.i1.AND.lbndg(1,j).EQ.i2) ok=.false.
              IF (lbndg(2,j).EQ.i1.AND.lbndg(3,j).EQ.i2) ok=.false.
              IF (lbndg(3,j).EQ.i1.AND.lbndg(2,j).EQ.i2) ok=.false.
          END DO
          do itors=1,n
             if((i1.eq.int14(1,itors).and. i2.eq.int14(2,itors)).or.
     &            (i2.eq.int14(1,itors).and. i1.eq.int14(2,itors))) THEN
                ok = .false.
             end if   
          end do
          IF(ok) THEN
              n=n+1
              int14(1,n)=i1
              int14(2,n)=i2
          END IF
          IF(.NOT. dihed) THEN
              itorsl(1,i)=0
              itorsl(2,i)=0
              itorsl(3,i)=0
              itorsl(4,i)=0
          END IF
      END DO

      int14p=n
      IF(.NOT. dihed) THEN
          ntors=0
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
