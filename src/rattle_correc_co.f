      SUBROUTINE rattle_correc_co(co,dss,cnst,vco,mass,degree,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <98/02/12 15:57:31 marchi>                             *
*                                                                      *
*   Constraint the co matrix to support only isotropic changes         *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
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
      REAL*8  co(3,3),vco(3,3),mass(*),dss(5)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,iox,iter
      REAL*8 tol
      REAL*8  dpp,gg,dcnst(2,5),cc(5)
      DATA tol/1.0D-9/

*==================== EXECUTABLE STATEMENTS ============================


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
*---- RATTLE loop ------------------------------------------------------
*=======================================================================

      iter=0
      iox=1
      DO WHILE(iox .NE. 0)
         iox=0
         DO k=1,degree
            la=cnst(1,k)
            lb=cnst(2,k)
            IF(DABS(dss(k)) .LT. 1.0D-6) THEN
               vco(la,lb)=0.0D0
            ELSE
               dpp=vco(la,lb)-dss(k)*vco(1,1)
               cc(k)=dpp
               IF(DABS(dpp) .GT. tol) THEN
                  iox=1
                  gg=dcnst(1,k)/mass(la)-dcnst(2,k)/mass(1)
                  gg=dpp/gg
                  vco(la,lb)=vco(la,lb)-gg*dcnst(1,k)/mass(la)
                  vco(1,1)=vco(1,1)-gg*dcnst(2,k)/mass(1)
               END IF
            END IF
         END DO
         IF(iox.ne.0) THEN
            iter=iter+1
            IF(iter.GT.5000)THEN
               iret=1
               errmsg=
     & ' While RATTLEing CO matrix: The iteration procedure did not'
     x              //' converge.'
               RETURN
            END IF
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
