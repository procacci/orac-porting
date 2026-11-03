      SUBROUTINE prtit(beta,res,tors,torsp,ptors,o1)

************************************************************************
*                                                                      *
*     PRTIT will print the improper torsion initial angles and         *
*     potential multiplicities.                                        *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER o1,torsp,tors(4,*),res(*)
      REAL*8  co(3,3),ptors(o1,*)
      CHARACTER*7 beta(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,l1,l2,l3,l4,ntph
      REAL*8  coa

*==================== EXECUTABLE STATEMENTS ============================

      IF(torsp .EQ. 0) RETURN
      WRITE(kprint,100)
      WRITE(kprint,200)
      DO i=1,torsp
         l1=tors(1,i)
         l2=tors(2,i)
         l3=tors(3,i)
         l4=tors(4,i)
         IF(ptors(i,3) .GT. 0.0D0) THEN
            coa=ptors(i,1)*efact/1000.0D0
            WRITE(kprint,300) i,beta(l1),res(l1),beta(l2),res(l2)
     &           ,beta(l3),res(l3),beta(l4),res(l4),coa,ptors(i,2)*180
     &           .0D0/pi
         ELSE
            coa=ptors(i,1)*efact/1000.0D0
            ntph=DINT(ptors(i,2)+0.5D0)
            WRITE(kprint,400) i,beta(l1),l1,beta(l2),l2
     &           ,beta(l3),l3,beta(l4),l4,coa,ntph
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(//1x,'*',75('='),'*'/1x,'*',10('='),
     x '    E q u i l b .   I - T o r s i o n   A n g l e s    ',
     x 10('='),'*'/1x,'*',75('='),'*'/)
200   FORMAT(/5x,'  Atom1( No  ) ','  Atom2( No  ) ','  Atom3( No  ) ',
     x'  Atom4( No  ) ',4x,' Force ',1x,' 0Angle/Mult. '/)
300   FORMAT(1x,i5,1x,a5,'(',i5,')',3x,a5,'(',i5,')',3x,a5,'(',i5,')',3x
     &     ,a5,'(',i5,')',1x,f10.2,f10.2)
400   FORMAT(1x,i5,1x,a5,'(',i5,')',3x,a5,'(',i5,')',3x,a5,'(',i5,')',3x
     &     ,a5,'(',i5,')',1x,f10.2,2x,i5)
      RETURN
      END
