      SUBROUTINE prtba(beta,res,bend,bendp,pbend,o1)

************************************************************************
*                                                                      *
*     PRTBA will print the bond angle initial and equilibrium          *
*     values.                                                          *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER o1,bendp,bend(3,*),res(*)
      REAL*8  pbend(o1,*)
      CHARACTER*7 beta(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,lc
      REAL*8  b0

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kprint,100)
      WRITE(kprint,200)
      DO i=1,bendp
         la=bend(1,i)
         lb=bend(2,i)
         lc=bend(3,i)
         b0=pbend(i,2)*180.0/pi
         WRITE(kprint,300) i,beta(la),la,beta(lb),lb,beta(lc),lc,pbend(i
     &        ,1)*efact/1000.0D0,b0
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(//1x,'*',75('='),'*'/1x,'*',15('='),
     x '    E q u i l .     B o n d   A n g l e s    ',15('='),'*'/
     x     1x,'*',75('='),'*'/)
200   FORMAT(/7x,'  Atom1( No  ) ','  Atom2( No  ) ','  Atom3( No  ) ',
     x 6x,' Force  ',6x,' 0Angle '/)
300   FORMAT(3x,i5,2x,a5,'(',i5,')',4x,a5,'(',i5,')',4x,a5,'(',i5,')',4x
     &     ,f10.2,2x,f10.2,f10.2)
      RETURN
      END
