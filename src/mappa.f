      SUBROUTINE mappa(bond,lbond,bend,lbend,int14,int14p,nato,map,
     x                na,lwork)

************************************************************************
*   Time-stamp: <97/09/12 12:43:43 marchi>                             *
*                                                                      *
*                                                                      *
*     Calculate the non bonded interaction map. Include in the map     *
*     the adddresses for atoms interacting in 1-2, 1-3 and 1-4.        *
*                                                                      *
*     BEND    :  List of bendings for the macromolecule.      (INPUT)  *
*                >> integer BEND(3,*) <<                               *
*     LBEND   :  Number of bendings.                          (INPUT)  *
*     INT14   :  List of all 1-4 interactions.                (INPUT)  *
*                >> integer INT14(2,*) <<                              *
*     INT14P  :  Number of the 1-4 interactions.              (INPUT)  *
*     NATO    :  Number of atoms forming the macromolecule.   (INPUT)  *
*     MAP     :  Non bonded interaction map.                 (OUTPUT)  *
*                >> integer   MAP(*) <<                                *
*                See START heading for explanations.                   *
*     LWORK   :  Logical work array.                          (INPUT)  *
*                >> logical*4 LWORK(*) <<                              *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jul  2 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n1,na,lbond,lbend,int14p,nato,bond(2,*),bend(3,*),
     x        int14(2,*)
      INTEGER map(*)
      LOGICAL lwork(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,j1,j2,l,n,m,noff
      LOGICAL ok

*==================== EXECUTABLE STATEMENTS ============================


      n=0
      DO i=1,nato
         DO j=i+1,nato
            lwork(j-i)=.TRUE.
         END DO
         DO l=1,lbond
            ok=.FALSE.
            IF(bond(1,l) .EQ. i) THEN
               j1=bond(2,l)
               ok=.TRUE.
            END IF
            IF(bond(2,l) .EQ. i) THEN
               j1=bond(1,l)
               ok=.TRUE.
            END IF
            IF(ok) THEN
               IF(j1 .GT. i) THEN
                  lwork(j1-i)=.FALSE.
               END IF
            END IF
         END DO
         

         DO l=1,lbend
            ok=.FALSE.
            IF(bend(1,l).EQ.i) THEN
               j1=bend(2,l)
               j2=bend(3,l)
               ok=.TRUE.
            END IF
            IF(bend(2,l).EQ.i) THEN
               j1=bend(3,l)
               j2=bend(1,l)
               ok=.TRUE.
            END IF
            IF(bend(3,l).EQ.i) THEN
               j1=bend(1,l)
               j2=bend(2,l)
               ok=.TRUE.
            END IF
            IF(ok) THEN
               IF(j1.GT.i) THEN
                  lwork(j1-i)=.FALSE.
               END IF
               IF(j2.GT.i) THEN
                  lwork(j2-i)=.FALSE.
               END IF
            END IF
         END DO
         
         IF(int14p.NE.0) THEN
            DO l=1,int14p
               IF(int14(1,l).EQ.i) THEN
                  j1=int14(2,l)
                  IF(j1.GT.i) THEN
                     lwork(j1-i)=.FALSE.
                  END IF
               ELSE IF(int14(2,l).EQ.i) THEN
                  j1=int14(1,l)
                  IF(j1.GT.i) THEN
                     lwork(j1-i)=.FALSE.
                  END IF
               END IF
            END DO
         END IF
         m=0
         DO j=i+1,nato
            IF(.NOT.lwork(j-i)) THEN
               m=m+1
               map(m+n+1)=j
            END IF
         END DO
         map(1+n)=m
         noff=m+1
         n=n+noff
         IF(na.LT.n) THEN
            errmsg=' IN MAPPA : Dimensions of MAPNL '//
     &           'are insufficient. ABORT !! '
            WRITE(6,'(10h NCOUNT = ,i10)') n
            CALL xerror(errmsg,80,1,2)
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
