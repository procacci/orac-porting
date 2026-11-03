      SUBROUTINE rattle_verlet_co(dt,co1,co,dss,cnst,vco,mass,degree
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <98/02/12 15:57:43 marchi>                             *
*                                                                      *
*   Constraint the co matrix to support only isotropic changes         *
*                                                                      *
*   The following constraints are satisfied:                           *
*                                                                      *
*        \sigma_{1} = {h_{22} \over h_{11}} - C_1                      *
*        \sigma_{2} = {h_{33} \over h_{11}} - C_2                      *
*        \sigma_{3} = {h_{12} \over h_{11}} - C_3                      *
*        \sigma_{4} = {h_{13} \over h_{11}} - C_4                      *
*        \sigma_{5} = {h_{23} \over h_{11}} - C_5                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 13 1997 -                                     *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,cnst(2,5),degree
      REAL*8  co1(3,3),co(3,3),vco(3,3),mass(*),dss(5),dt
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,iox,iter
      REAL*8 tol
      REAL*8  dcnt,dti,dp1,dp2,dpp,gg,dps1,dps2,dcnst(2,5)
      DATA tol/1.0D-9/

*==================== EXECUTABLE STATEMENTS ============================


      dti=1.0D0/dt

*=======================================================================
*---- Compute the vectors associated with each constraint and store ----
*---- them in a temporary array ----------------------------------------
*=======================================================================

      DO k=1,degree
         la=cnst(1,k)
         lb=cnst(2,k)
         dcnst(1,k)=1.0D0/co(1,1)
         dcnst(2,k)=-co(la,lb)/co(1,1)**2
      END DO

*=======================================================================
*---- SHAKE loop -------------------------------------------------------
*=======================================================================

      iter=0
      iox=1
      DO WHILE(iox .NE. 0)
         iox=0
         DO k=1,degree
            la=cnst(1,k)
            lb=cnst(2,k)
            IF(DABS(dss(k)) .LT. 1.0D-6) THEN
               co1(la,lb)=co(la,lb)
            ELSE
               dp1=1.0D0/co1(1,1)
               dp2=-co1(la,lb)/co1(1,1)
               dpp=co1(la,lb)/co1(1,1)
               dcnt=dpp-dss(k)
               IF(DABS(dcnt) .GT. tol) THEN
                  iox=1
                  dps1=dp1*dcnst(1,k)
                  dps2=dp2*dcnst(2,k)
                  gg=dps1/mass(la)+dps2/mass(1)
                  gg=-dcnt/gg
                  co1(la,lb)=co1(la,lb)+dcnst(1,k)*gg/mass(la)
                  co1(1,1)=co1(1,1)+dcnst(2,k)*gg/mass(1)

                  gg=gg*dti

                  vco(la,lb)=vco(la,lb)+dcnst(1,k)*gg/mass(la)
                  vco(1,1)=vco(1,1)+dcnst(2,k)*gg/mass(1)
               END IF
            END IF
         END DO
         IF(iox.ne.0) THEN
            iter=iter+1
            IF(iter.GT.5000)THEN
               iret=1
               errmsg=
     & ' While SHAKEing CO matrix: The iteration procedure did not'
     x              //' converge.'
               RETURN
            END IF
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
