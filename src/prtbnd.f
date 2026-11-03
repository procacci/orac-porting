      SUBROUTINE prtbnd(beta,res,bond,bondp,pbond,o1)

************************************************************************
*                                                                      *
*     PRTBA will print the bond stretching                             *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER o1,bondp,bond(2,*),res(*)
      REAL*8  pbond(o1,*)
      CHARACTER*7 beta(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,lc
      REAL*8 b0

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kprint,100)
      WRITE(kprint,200)
      DO i=1,bondp
         la=bond(1,i)
         lb=bond(2,i)
         b0=pbond(i,2)
         WRITE(kprint,300) i,beta(la),la,beta(lb),lb,pbond(i,1)*efact
     &        /1000.0D0,b0
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(//1x,'*',75('='),'*'/1x,'*',15('='),
     x '    B o n d   S t r e t c h i n g s          ',15('='),'*'/
     x     1x,'*',75('='),'*'/)
200   FORMAT(/7x,'  Atom1( No  ) ','  Atom2( No  ) ',
     x 6x,' Force ',6x,' 0Bond '/)
300   FORMAT(3x,i5,2x,a5,'(',i5,')',4x,a5,'(',i5,')',4x,f10.2,f10.4)
      RETURN
      END
