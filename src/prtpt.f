      SUBROUTINE prtpt(beta,res,tors,torsp,ptors,o1)

************************************************************************
*                                                                      *
*     PRTPT will print the proper torsions initial angles and the      *
*     potential multiplicities.                                        *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER torsp,o1,tors(4,*),res(*)
      CHARACTER*7 beta(*)
      REAL*8  ptors(o1,*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,l1,l2,l3,l4,ntph
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,coa,aux

*==================== EXECUTABLE STATEMENTS ============================

      IF(torsp .EQ. 0) RETURN
      WRITE(kprint,100)
      WRITE(kprint,200)
      DO 10 i=1,torsp
          l1=tors(1,i)
          l2=tors(2,i)
          l3=tors(3,i)
          l4=tors(4,i)
          coa=ptors(i,1)*efact/1000.0D0
          ntph=DINT(ptors(i,2)+0.5D0)
          WRITE(kprint,300) i,beta(l1),l1,beta(l2),l2,beta(l3)
     &         ,l3,beta(l4),l4,coa,ntph
10    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(//1x,'*',75('='),'*'/1x,'*',10('='),
     x '    I n i t i a l   P - T o r s i o n   A n g l e s    ',
     x 10('='),'*'/1x,'*',75('='),'*'/)
200   FORMAT(/5x,'  Atom1( No  ) ','  Atom2( No  ) ','  Atom3( No  ) ',
     x'  Atom4( No  ) ',4x,' Force ',1x,' Mult '/)
300   FORMAT(1x,i5,1x,a5,'(',i5,')',3x,a5,'(',i5,')',3x,a5,'(',i5,')',3x
     &     ,a5,'(',i5,')',1x,f10.2,1x,i5)
      RETURN
      END
