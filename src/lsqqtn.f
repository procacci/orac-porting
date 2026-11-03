      SUBROUTINE lsqqtn(n,r,r0,w,q,rfit,error,iopt,iret,msg)
 
**********************************************************************
*                                                                    *
*     SR LSQQTN fits a rigid structure of N sites into a given       *
*     structure of N sites by a ROTATION. All coordinates are        *
*     related the origin/center of mass ZERO.                        *
*                                                                    *
*     AUTHOR: Gerald Kneller                                         *
*             IBM France                                             *
*             94-96, rue Reaumur                                     *
*             F-75002 Paris, FRANCE                                  *
*                                                                    *
*---- LAST UPDATE 20.07.90 ------------------------------------------*
*                                                                    *
*     ARGUMENTS:                                                     *
*                                                                    *
*     N      : Number of sites.                                 (I)  *
*              >> INTEGER N <<                                       *
*     R      : Input coordinates.                               (I)  *
*              >> REAL*8 R(3,*) <<                                   *
*     R0     : Definition of the rigid structure.               (I)  *
*              >> REAL*8 R0(3,*) <<                                  *
*     W      : Relative weights of the sites.                   (I)  *
*              >> REAL*8 W(*) <<                                     *
*     Q      : Fitted quaternion parameters.                    (O)  *
*              If the solution is unique, they are are stored        *
*              in q(0),..,q(3).                                      *
*              If the solution is not unique, the two linear         *
*              quaternion-vectors are stored in q(0),..,q(3) and     *
*              q(4),..,q(7).                                         *
*              >> REAL*8 Q(0:7) <<                                   *
*     RFIT   : Fitted coordinates  (IOPT <> 0).                 (O)  *
*              >> REAL*8 RFIT(3,*) <<                                *
*     ERROR  : RMS error / site.                                (O)  *
*              >> REAL*8 ERROR <<                                    *
*     IOPT   : Options:                                         (I)  *
*              IOPT =  0: only quaternions are calculated,           *
*              IOPT <> 1: quaternions and fitted coodinates          *
*                         are calculated. RFIT and R can be          *
*                         equivalenced, in this case R will be       *
*                         overwritten by RFIT on output.             *
*              >> INTEGER IOPT <<                                    *
*     IRET   : Return-code:                                     (O)  *
*              IRET =  0: normal execution, unique solution.         *
*              IRET = 10: normal execution, non-unique solution.     *
*              IRET =  9: bad performance of EIGRS, unique solution. *
*              IRET = 19: bad perform. of EIGRS, non-unique solution.*
*              IRET = 99: Fatal error in EIGRS.                      *
*              >> INTEGER iret <<                                    *
*     MSG    : Message corresponding to IRET.                   (O)  *
*              >> CHARACTER*80 msg <<                                *
*                                                                    *
*     EXTERNALS:                                                     *
*                                                                    *
*     - subroutine EIGRS (IMSL)                                      *
*                                                                    *
**********************************************************************
 
*==== DECLARATIONS: =================================================*
 
      IMPLICIT NONE
 
*---- ARGUMENTS: ----------------------------------------------------*
 
      integer      n,iopt,iret
      real*8       r(3,*),r0(3,*),rfit(3,*),w(*),error,q(0:7)
      character*80 msg
 
*---- LOCAL VARIABLES: ----------------------------------------------*
 
      integer i,j,jobn,ier
      real*8  epsi,s,
     .        m(4,4),d(3,3),lambda(4),z(4,4),wk(20)
 
*---- DATA STATEMENTS: ----------------------------------------------*
 
      data epsi /1.0e-10/
 
*==== EXECUTABLE STATEMENTS: ========================================*
 
      iret = 0
 
*---- MATRIX M FOR EIGENVALUE PROBLEM: ------------------------------*
 
      do 100 i = 1,4
         do 100 j = 1,4
            m(i,j) = 0.d0
100   continue
 
      do 200 i = 1,n
         if (w(i) .eq. 0.d0) goto 200
         m(1,1)=m(1,1)
     .    +w(i)*(  r0(1,i)*r0(1,i) + r0(2,i)*r0(2,i) + r0(3,i)*r0(3,i)
     .           + r(1,i)*r(1,i)   + r(2,i)*r(2,i)   +  r(3,i)*r(3,i)
     .           - 2.d0*r0(1,i)*r(1,i)
     .           - 2.d0*r0(2,i)*r(2,i)
     .           - 2.d0*r0(3,i)*r(3,i) )
         m(1,2)=m(1,2)
     .             +2.d0*w(i)*(-r0(2,i)*r(3,i)+r0(3,i)*r(2,i))
         m(1,3)=m(1,3)
     .             +2.d0*w(i)*(r0(1,i)*r(3,i)-r0(3,i)*r(1,i))
         m(1,4)=m(1,4)
     .             +2.d0*w(i)*(-r0(1,i)*r(2,i)+r0(2,i)*r(1,i))
         m(2,2)=m(2,2)
     .    +w(i)*(  r0(1,i)*r0(1,i) + r0(2,i)*r0(2,i) + r0(3,i)*r0(3,i)
     .           + r(1,i)*r(1,i)   + r(2,i)*r(2,i)   +  r(3,i)*r(3,i)
     .           - 2.d0*r0(1,i)*r(1,i)
     .           + 2.d0*r0(2,i)*r(2,i)
     .                    + 2.d0*r0(3,i)*r(3,i) )
         m(2,3)=m(2,3)
     .             -2.d0*w(i)*(r0(1,i)*r(2,i)+r0(2,i)*r(1,i))
         m(2,4)=m(2,4)
     .          -2.d0*w(i)*(r0(1,i)*r(3,i)+r0(3,i)*r(1,i))
         m(3,3)=m(3,3)
     .    +w(i)*(  r0(1,i)*r0(1,i) + r0(2,i)*r0(2,i) + r0(3,i)*r0(3,i)
     .           + r(1,i)*r(1,i)   + r(2,i)*r(2,i)   +  r(3,i)*r(3,i)
     .           + 2.d0*r0(1,i)*r(1,i)
     .           - 2.d0*r0(2,i)*r(2,i)
     .           + 2.d0*r0(3,i)*r(3,i) )
         m(3,4)=m(3,4)
     .             -2.d0*w(i)*(r0(2,i)*r(3,i)+r0(3,i)*r(2,i))
         m(4,4)=m(4,4)
     .    +w(i)*(  r0(1,i)*r0(1,i) + r0(2,i)*r0(2,i) + r0(3,i)*r0(3,i)
     .           + r(1,i)*r(1,i)   + r(2,i)*r(2,i)   +  r(3,i)*r(3,i)
     .           + 2.d0*r0(1,i)*r(1,i)
     .           + 2.d0*r0(2,i)*r(2,i)
     .           - 2.d0*r0(3,i)*r(3,i) )
200   continue
 
      m(2,1) = m(1,2)
 
      m(3,1) = m(1,3)
      m(3,2) = m(2,3)
 
      m(4,1) = m(1,4)
      m(4,2) = m(2,4)
      m(4,3) = m(3,4)
 
*---- SOLVE THE EIGENVECTOR PROBLEM FOR M: --------------------------*
 
      jobn   = 12
      call EIGRS(m,4,jobn,lambda,z,4,wk,ier)
      if (ier .ne. 0)      goto 9000
      if (wk(1) .gt. 1.d0) iret = 9
 
*---- PICK THE CORRECT EIGENVECTOR(S): ------------------------------*
 
      if (z(1,1) .lt. 0.d0) then
         s = -1.d0
      else
         s =  1.d0
      endif
      q(0) = s*z(1,1)
      q(1) = s*z(2,1)
      q(2) = s*z(3,1)
      q(3) = s*z(4,1)
      if (dabs(lambda(1)) .lt. epsi) then
         error = 0.d0
      else
         error = dsqrt(lambda(1))
      endif
 
      if (dabs(lambda(1) - lambda(2)) .lt. epsi) then
         iret = iret + 10
         if (z(1,2) .lt. 0.d0) then
            s = -1.d0
         else
            s =  1.d0
         endif
         q(4) = s*z(1,2)
         q(5) = s*z(2,2)
         q(6) = s*z(3,2)
         q(7) = s*z(4,2)
      else
         q(4) = 0.d0
         q(5) = 0.d0
         q(6) = 0.d0
         q(7) = 0.d0
      endif
 
*.....RETURN IF NO FITTED COORDINATES ARE NEEDED:
      if (iopt .eq. 0) goto 8000
 
*---- ROTATION MATRIX IN TERMS OF QUATERNIONS: ----------------------*
 
      d(1,1)=-2.d0*q(2)**2-2.d0*q(3)**2+1.d0
      d(1,2)=2.d0*(-q(0)*q(3)+q(1)*q(2))
      d(1,3)=2.d0*(q(0)*q(2)+q(1)*q(3))
      d(2,1)=2.d0*(q(0)*q(3)+q(1)*q(2))
      d(2,2)=-2.d0*q(1)**2-2.d0*q(3)**2+1.d0
      d(2,3)=2.d0*(-q(0)*q(1)+q(2)*q(3))
      d(3,1)=2.d0*(-q(0)*q(2)+q(1)*q(3))
      d(3,2)=2.d0*(q(0)*q(1)+q(2)*q(3))
      d(3,3)=-2.d0*q(1)**2-2.d0*q(2)**2+1.d0
 
*---- CALCULATE FIT: ------------------------------------------------*
 
      do 400 i = 1,n
         rfit(1,i) = d(1,1)*r0(1,i) + d(1,2)*r0(2,i) + d(1,3)*r0(3,i)
         rfit(2,i) = d(2,1)*r0(1,i) + d(2,2)*r0(2,i) + d(2,3)*r0(3,i)
         rfit(3,i) = d(3,1)*r0(1,i) + d(3,2)*r0(2,i) + d(3,3)*r0(3,i)
400   continue
 
*---- RETURN: -------------------------------------------------------*
 
8000  if (iret .eq. 0) then
         msg = ' LSQQTN: Normal execution, unique solution'
         return
      endif
 
      if (iret .eq. 10) then
         msg = ' LSQQTN: Normal execution, non-unique solution'
         return
      endif
 
      if (iret .eq. 9) then
         msg = ' LSQQTN: Bad perform. in eigrs, unique solution'
         return
      endif
 
      if (iret .eq. 19) then
         msg = ' LSQQTN: Bad perform. in eigrs, non-unique solution'
         return
      endif
 
9000  iret = 99
      msg = ' LSQQTN: Fatal error in eigrs'
 
      return
      end
