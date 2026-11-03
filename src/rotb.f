      SUBROUTINE rotb(a,b,c,alfa,beta,gamma,co)

************************************************************************
*                                                                      *
*                                                                      *
*     A       :  Unit cell length along the x direction.               *
*     B       :  Unit cell length along the y direction.               *
*     C       :  Unit cell length along the z direction.               *
*     ALFA    :  Angle between the y and z axis.                       *
*     BETA    :  Angle between the z and x axis.                       *
*     GAMMA   :  Angle between the x and y axis.                       *
*     CO      :  Rotation matrix from simulation box frame to          *
*                orthogonalized frame.                                 *
*                >> real*8 CO(3,3) <<                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 a,b,c,alfa,beta,gamma,co(3,3)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  degrad,qt,alf,bet,gam,ax,bx,by,cx,cy,cz

*==================== EXECUTABLE STATEMENTS ============================

      degrad=pi/180.0d0

      ax=co(1,1)*boxl
      bx=co(1,2)*boxl
      cx=co(1,3)*boxl
      by=co(2,2)*boxl
      cy=co(2,3)*boxl
      cz=co(3,3)*boxl

      a=ax
      b=DSQRT(bx**2+by**2)
      c=DSQRT(cx**2+cy**2+cz**2)
      gam=bx/b
      qt=by/b
      bet=cx/c
      alf=cy*qt/c+bet*gam

      alfa=DACOS(alf)/degrad
      beta=DACOS(bet)/degrad
      gamma=DACOS(gam)/degrad
      qt=DASIN(qt)/degrad

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
