      SUBROUTINE addh2o(x0,y0,z0)

************************************************************************
*                                                                      *
*     ADDH2O will produce the coordinates of the hydrogens for         *
*     a water molecule.                                                *
*                                                                      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(2), Y0(2), Z0(2) <<                      *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL COMD0                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0(*),y0(*),z0(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(2),y(2),z(2),vx,vy,vz,ux,uy,uz,ranf
      REAL*8 norm,pi,angle,length,alpha,eps,beta
      INTEGER m
      external ranf

!=======================================================================
*---- Parameters for SPC water -----------------------------------------
*=======================================================================

      DATA angle, length, eps/109.47D0, 1.0D0, 1.0D-8/

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute a first randomly oriented versor -------------------------
*=======================================================================

      x(1)=2.0D0*(ranf()-0.5D0)
      y(1)=2.0D0*(ranf()-0.5D0)
      z(1)=2.0D0*(ranf()-0.5D0)
      norm=x(1)**2+y(1)**2+z(1)**2
      norm=DSQRT(norm)
      x(1)=x(1)/norm
      y(1)=y(1)/norm
      z(1)=z(1)/norm

*=======================================================================
*---- Compute the second versor at angle "angle" from the first --------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      alpha=angle*pi/180.0D0
      beta=(180-angle)*pi/180.0D0
      vx=2.0D0*(ranf()-0.5D0)
      vy=2.0D0*(ranf()-0.5D0)
      alpha=DCOS(alpha)
      vz=-vx*x(1)-vy*y(1)
      vz=vz/z(1)
      norm=vx**2+vy**2+vz**2
      norm=DSQRT(norm)
      vx=DSIN(beta)*vx/norm
      vy=DSIN(beta)*vy/norm
      vz=DSIN(beta)*vz/norm
      ux=-x(1)*DCOS(beta)
      uy=-y(1)*DCOS(beta)
      uz=-z(1)*DCOS(beta)
      x(2)=vx+ux
      y(2)=vy+uy
      z(2)=vz+uz
      norm=x(2)**2+y(2)**2+z(2)**2
      norm=DSQRT(norm)
      x(2)=x(2)/norm
      y(2)=y(2)/norm
      z(2)=z(2)/norm

*=======================================================================
*---- Create the two new vector ----------------------------------------
*=======================================================================

      DO m=1,2
          x0(m)=x(m)*length
          y0(m)=y(m)*length
          z0(m)=z(m)*length
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
