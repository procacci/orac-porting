      SUBROUTINE bcnstp(ss_index,cnst,ncnst,xcoor,ycoor,zcoor,mass
     &     ,dss,coeff,cnstp,ma,cnstpp,cnstpp_slv,iret,errmsg)

************************************************************************
*                                                                      *
*     Compute two vectors which are used by the constraints            *
*     subroutine.                                                      *
*                                                                      *
*                                                                      *
*     CNST    :  Constraint matrix.                           (INPUT)  *
*                >> integer CNST(2,N1) <<                              *
*     NCNST   :  Number of constraints.                       (INPUT)  *
*     N1      :  Physical dimension of the constraint                  *
*                matrix.                                               *
*     XCOOR   :  Atomic cartesian coordinates.                (INPUT)  *
*                >> real*8 XCOOR(*) <<                                 *
*     YCOOR   :  Atomic cartesian coordinates.                (INPUT)  *
*                >> real*8 YCOOR(*) <<                                 *
*     ZCOOR   :  Atomic cartesian coordinates.                (INPUT)  *
*                >> real*8 ZCOOR(*) <<                                 *
*     NATOM   :  Number of atoms involvend in the             (INPUT)  *
*                constraint.                                           *
*     MASS    :  Atomic mass.                                 (INPUT)  *
*     DSS     :  Square od the distance associated with      (OUTPUT)  *
*                each constraint.                                      *
*                >> real*8 DSS(*) <<                                   *
*     COEFF   :  Coefficient used by the constraint          (OUTPUT)  *
*                algorithm.                                            *
*                >> real*8 COEFF(*) <<                                 *
*                                                                      *
*---- Last update 07/08/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,ma,ncnst,cnstpp,cnstpp_slv,cnst(2,*),cnstp(2,*)
     &     ,ss_index(*)
      REAL*8  xcoor(*),ycoor(*),zcoor(*),mass(*)
      REAL*8  dss(*),coeff(*)
      CHARACTER*80 errmsg
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER k,la,lb,noff,m
      REAL*8 dsx,dsy,dsz

*==================== EXECUTABLE STATEMENTS ============================

      noff=cnstpp
      IF(noff+ncnst .GT. ma) THEN
          iret=1
          errmsg='In BCNSTP: Number of constraints exceeds physical'
     x//         ' dimension. Abort.'
          RETURN
      END IF
      cnstpp_slv=0 
      DO k=1,ncnst
          noff=noff+1
          la=cnst(1,k)
          lb=cnst(2,k)
          m=ss_index(la)
          dsx=xcoor(la)-xcoor(lb)
          dsy=ycoor(la)-ycoor(lb)
          dsz=zcoor(la)-zcoor(lb)
          dss(noff)=dsx**2+dsy**2+dsz**2
          IF(m .EQ. 2) cnstpp_slv=cnstpp_slv+1
          IF(.NOT. (near0(mass(la)) .OR. near0(mass(lb)))) THEN
              coeff(noff)=2.0D0*(mass(la)+mass(lb))/
     x              (mass(la)*mass(lb))
          ELSE
              coeff(noff)=0.0D0
          END IF
          cnstp(1,noff)=la
          cnstp(2,noff)=lb
      END DO
      cnstpp=noff

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
