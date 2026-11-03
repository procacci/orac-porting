      SUBROUTINE comd0(a,b,c,alpha,d0)

************************************************************************
*                                                                      *
*                                                                      *
*     COMD0 will calculate the rotation matrix D0 which rotate         *
*     a vector onto the direction CHI=(a,b,c) and apply a clockwise    *
*     rotation of alpha radians.                                       *
*                                                                      *
*     A       :  Components of the reference versor.              (I)  *
*     B                                                                *
*     C                                                                *
*     ALPHA   :  Clockwise rotation angle about (a,b,c).          (I)  *
*     D0      :  Rotation matrix.                                 (O)  *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT CHARACTER*80(a-z)

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 a,b,c,alpha,d0(3,3)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 d,cs,sn

*==================== EXECUTABLE STATEMENTS ============================

      d=DSQRT(a**2+b**2)
      IF(DABS(d).GE.1.0D-8) THEN
          cs=DCOS(alpha)
          sn=DSIN(alpha)
          d0(1,1)=(b*cs-a*c*sn)/d
          d0(1,2)=(b*sn+a*c*cs)/d
          d0(1,3)=a
          d0(2,1)=-(a*cs+b*c*sn)/d
          d0(2,2)=(-a*sn+b*c*cs)/d
          d0(2,3)=b
          d0(3,1)=sn*d
          d0(3,2)=-cs*d
          d0(3,3)=c
      ELSE
          cs=DCOS(alpha)
          sn=DSIN(alpha)
          d0(1,1)=cs
          d0(1,2)=sn
          d0(1,3)=0.0D0
          d0(2,1)=-sn
          d0(2,2)=cs
          d0(2,3)=0.0D0
          d0(3,1)=0.0D0
          d0(3,2)=0.0D0
          d0(3,3)=1
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
