      SUBROUTINE covbd(potbo,potbe,potto,potit,ptorj,lbond,lbend,ltors,
     x                 litor,itor_ptype,m0,m1,m2,m3)

************************************************************************
*                                                                      *
*     This subroutine converts the covalent interaction parameters     *
*     associated to the lists of bendings, torsions and improper       *
*     torsions to program units.                                       *
*                                                                      *
*     POTBE   :  List of bending potentials.          (INPUT/OUTPUT)   *
*                >> real*8 POTBE(M1,2) <<                              *
*     POTTO   :  List of torsion potentials.          (INPUT/OUTPUT)   *
*                >> real*8 POTTO(M2,2) <<                              *
*     POTIT   :  List of improper torsion potentials. (INPUT/OUTPUT)   *
*                >> real*8 POTIT(M3,2) <<                              *
*     LBEND   :  Number of bendings.                          (INPUT)  *
*     LTORS   :  Number of proper torsions.                   (INPUT)  *
*     LITOR   :  Number of improper torsions.                 (INPUT)  *
*     M1      :  Physical dimension of POTBE.                 (INPUT)  *
*     M2      :  Physical dimension of POTTO.                 (INPUT)  *
*     M3      :  Physical dimension of POTIT.                 (INPUT)  *
*                                                                      *
*------ Last update 04/11/89 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER lbond,lbend,ltors,litor,m0,m1,m2,m3,itor_ptype
      REAL*8 potbo(m0,*),potbe(m1,*),potto(m2,*),potit(m3,*),ptorj(m2,*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j

*==================== EXECUTABLE STATEMENTS ============================

      DO 5 i=1,lbond
          potbo(i,1)=1000.0D0*potbo(i,1)*4.184/(unite*avogad)
5     CONTINUE
      DO 10 i=1,lbend
          potbe(i,1)=1000.0D0*potbe(i,1)*4.184/(unite*avogad)
          potbe(i,2)=potbe(i,2)*(pi/180.0D0)
          potbe(i,3)=1000.0D0*potbe(i,3)*4.184/(unite*avogad)
c--       if 7 bending parameters are different from zero then 
c--       it means  ka  a0  kub rub k12 r1  r2 
          if(abs(potbe(i,6)).gt.0.d0) THEN 
            potbe(i,5)=1000.0D0*potbe(i,5)*4.184/(unite*avogad)
          endif
10    CONTINUE
      DO 20 i=1,ltors
          potto(i,1)=1000.0D0*potto(i,1)*4.184/(unite*avogad)
          DO 40 j=1,4
              ptorj(i,j)=1000.0D0*ptorj(i,j)*4.184/(unite*avogad)
40        CONTINUE
20    CONTINUE
      IF(itor_ptype .NE. 0) THEN
         IF(itor_ptype .EQ. 1) THEN
            DO i=1,litor
               potit(i,1)=1000.0D0*potit(i,1)*4.184/(unite*avogad)
               potit(i,2)=potit(i,2)*(pi/180.0D0)
            END DO
         ELSE IF(itor_ptype .EQ. 2) THEN
            DO i=1,litor
               potit(i,1)=1000.0D0*potit(i,1)*4.184/(unite*avogad)
            END DO
         END IF
      ELSE
         DO i=1,litor
            potit(i,1)=1000.0D0*potit(i,1)*4.184/(unite*avogad)
            IF(potit(i,3) .GT. 0.0D0) THEN
               potit(i,2)=potit(i,2)*(pi/180.0D0)
            END IF
         END DO
      END IF
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
