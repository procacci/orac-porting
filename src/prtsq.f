      SUBROUTINE prtsq(nbun,mend,prsymb)

************************************************************************
*   Time-stamp: <97/07/14 08:17:36 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jul 13 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nbun,mend(*)
      CHARACTER*8 prsymb(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER index(n5),index2(nores+1,n5)
      COMMON /rag1/ index,index2

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,i1,n,map,ia,ib,j,m
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(kprint,300)
      DO i=1,nbun,7
         n=MIN(6,nbun-i)
         WRITE(kprint,100) (prsymb(mend(i1)),i1=i,i+n)
         WRITE(kprint,200) (i1,i1=i,i+n)
         WRITE(kprint,*)
      END DO

      WRITE(kprint,400)
      map=0
      DO i=1,nbun
         ia=mend(i)
         ok=.TRUE.
         DO j=1,map
            IF(ia .EQ. index(j)) ok=.FALSE.
         END DO
         IF(ok) THEN
            map=map+1
            index(map)=ia
         END IF
      END DO

      DO i=1,map
         index2(1,i)=0
         ia=index(i)
         DO j=1,nbun
            ib=mend(j)
            IF(ib .EQ. ia) THEN
               index2(1,i)=index2(1,i)+1
               m=index2(1,i)
               index2(1+m,i)=j
            END IF
         END DO
      END DO

      DO i=1,map
         m=index2(1,i)
         WRITE(kprint,500) prsymb(index(i)),m
         DO j=1,m,10
            n=MIN(9,m-j)
            WRITE(kprint,600) (index2(1+ia,i),ia=j,j+n)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(4x,7(a8,1x))
200   FORMAT(2x,7(i5,4x))
300   FORMAT(//1x,'*',75('='),'*'/1x,'*',16('='),
     x '       U n i t s   S e q u e n c e         ',16('='),'*'/
     x     1x,'*',75('='),'*'/)
400   FORMAT(//1x,'*',75('='),'*'/1x,'*',16('='),
     x '       U n i t s   p e r   T y p e         ',16('='),'*'/
     x     1x,'*',75('='),'*'/)
500   FORMAT(/10x,' Unit of type ',a8,' occurs ',i6,' times')
600   FORMAT(2x,10(i6,1x))

      RETURN
      END
