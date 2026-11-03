      SUBROUTINE add3h(ax,length,x0,y0,z0)

************************************************************************
*                                                                      *
*     ADD3H will produce coordinates for 3 hydrogens relative          *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(3), Y0(3), Z0(3) <<                      *
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

      REAL*8 x0(*),y0(*),z0(*),ax(3),length

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(3),y(3),z(3),d0(3,3)
      REAL*8 norm,alpha,beta,a,b,c,pi
      INTEGER i,m

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute the normalized versors -----------------------------------
*=======================================================================

      norm=0.0D0
      DO 10 i=1,3
          norm=norm+ax(i)**2
10    CONTINUE
      norm=DSQRT(norm)
      DO 20 i=1,3
          ax(i)=ax(i)/norm
20    CONTINUE

*=======================================================================
*---- Calculate the primitive coordinates of the hydrogens -------------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      beta=60.0D0*pi/180.0D0
      alpha=DACOS(-1.0D0/3.0D0)-0.50D0*pi
      x(1)=0.0D0
      x(2)= - DCOS(alpha)*DSIN(beta)
      x(3)=DCOS(alpha)*DSIN(beta)
      y(1)=-DCOS(alpha)
      y(2)=DCOS(alpha)*DCOS(beta)
      y(3)=y(2)
      z(1)=DSIN(alpha)
      z(2)=z(1)
      z(3)=z(1)
      a=ax(1)
      b=ax(2)
      c=ax(3)
      alpha=0.0D0
      CALL comd0(a,b,c,alpha,d0)

*=======================================================================
*---- Finally rotate the primitive coordinates and find the ------------
*---- required hydrogen coordinates ------------------------------------
*=======================================================================

      DO 30 m=1,3
          x0(m)=d0(1,1)*x(m)+d0(1,2)*y(m)+d0(1,3)*z(m)
          y0(m)=d0(2,1)*x(m)+d0(2,2)*y(m)+d0(2,3)*z(m)
          z0(m)=d0(3,1)*x(m)+d0(3,2)*y(m)+d0(3,3)*z(m)
30    CONTINUE
      DO 40 m=1,3
          x0(m)=x0(m)*length
          y0(m)=y0(m)*length
          z0(m)=z0(m)*length
40    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
