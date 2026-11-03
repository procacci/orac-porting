      SUBROUTINE rmbend(lbend,nl,mbend,nm)

************************************************************************
*                                                                      *
*     Remove angle bends from the list.                                *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*---- Last update 02/01/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley, CA 1992                   *
*                                                                      *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      use parst

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      INTEGER nl,nm
      INTEGER lbend(3,*),mbend(3,*)

!------------------ DEFINITION OF A SCRATCH COMMON ---------------------

      INTEGER bendx(3,m2)
      COMMON /rag1/ bendx

!-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,o,p,q,n,m,ncount
      LOGICAL ok

!==================== EXECUTABLE STATEMENTS ============================

      ncount=0
      DO m=1,nl
          o=lbend(1,m)
          p=lbend(2,m)
          q=lbend(3,m)
          DO n=1,nm
              i=mbend(1,n)
              j=mbend(2,n)
              k=mbend(3,n)
              ok=.NOT.((o .EQ. i .AND. p .EQ. j .AND. q .EQ. k) .OR.
     x           (o .EQ. k .AND. p .EQ. j .AND. q .EQ. i))
              IF(.NOT. ok) GOTO 100
          END DO
          ncount=ncount+1
          bendx(1,ncount)=lbend(1,m)
          bendx(2,ncount)=lbend(2,m)
          bendx(3,ncount)=lbend(3,m)
100       CONTINUE
      END DO
      nl=ncount
      DO n=1,nl
          lbend(1,n)=bendx(1,n)
          lbend(2,n)=bendx(2,n)
          lbend(3,n)=bendx(3,n)
      END DO

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
