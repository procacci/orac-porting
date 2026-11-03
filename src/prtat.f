      SUBROUTINE prtat(ss_index,nato,beta,betab,res,o1,prsymb,chrge,mass
     &     )

************************************************************************
*                                                                      *
*                                                                      *
*     Print atomic coordinates, types and charges.                     *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER o1
      INTEGER nato,res(o1,*),ss_index(*)
      REAL*8  chrge(*),mass(*)
      CHARACTER*8 prsymb(*)
      CHARACTER*7 beta(*),betab(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,typei
      CHARACTER*9  type

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kprint,100)
      WRITE(kprint,200)
      DO i=1,nato
         typei=ss_index(i)
         IF(typei .EQ. 1) type=' Solute  '
         IF(typei .EQ. 2) type=' Solvent '
         WRITE(kprint,300) i,beta(i),betab(i),res(i,1),prsymb(res(i,2))
     &        ,chrge(i)*DSQRT(unitc),mass(i),type
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(//1x,'*',75('='),'*'/1x,'*',16('='),
     x '    A t o m s   i n   t h e   S y s t e m  ',16('='),'*'/
     x     1x,'*',75('='),'*'/)
200   FORMAT('    No.  Label  Type        Residue            ',
     & 'Charge     Mass    Slt/Slv '/)
300   FORMAT(1x,i6,4x,a7,a7,i6,2x,a8,2x,f10.4,f10.4,2x,a9)

      RETURN
      END
