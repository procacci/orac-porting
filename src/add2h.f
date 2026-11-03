      SUBROUTINE add2h(ax,pax,iz,length,x0,y0,z0)

************************************************************************
*                                                                      *
*     ADD2H will produce coordinates for 2 hydrogens relative          *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     PAX     :  Axis perpendicular to the plane containing   (I)      *
*                the nearest neighbours of the central atom.           *
*                >> real*8 PAX(3) <<                                   *
*     IZ      :  It is 1 if PAX is defined and 0 otherwise.   (I)      *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(2), Y0(2), Z0(2) <<                      *
*                                                                      *
*---- Last update 07/08/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992.              *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0(*),y0(*),z0(*),ax(3),pax(3),length
      INTEGER iz

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(2),y(2),z(2),ver(3),vrr(3),vsg(3),d0(3,3)
      REAL*8 norm,alpha,a,b,c,pi,sgns
      INTEGER i,j,m

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute the normalized versors -----------------------------------
*=======================================================================

      norm=0.0D0
      DO i=1,3
          norm=norm+ax(i)**2
      END DO
      norm=DSQRT(norm)
      DO i=1,3
          ax(i)=ax(i)/norm
      END DO
      IF(iz.EQ.1) THEN
          norm=0.0D0
          DO i=1,3
              norm=norm+pax(i)**2
          END DO
          norm=DSQRT(norm)
          DO i=1,3
              pax(i)=pax(i)/norm
          END DO
      END IF

*=======================================================================
*---- Calculate the primitive coordinates of the hydrogens -------------
*---- and the versor perpendicular to their plane if required ----------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      alpha=0.50D0*(pi-DACOS(-1.0D0/3.0D0))
      IF(iz.NE.1) alpha=30.0D0*pi/180.0D0
      x(1)=0.0D0
      x(2)=0.0D0
      y(1)=DCOS(alpha)
      y(2)=-y(1)
      z(1)=DSIN(alpha)
      z(2)=z(1)
      a=ax(1)
      b=ax(2)
      c=ax(3)
      alpha=0.0D0
      CALL comd0(a,b,c,alpha,d0)
      IF(iz.EQ.1) THEN
          ver(1)=y(2)*z(1)-z(2)*y(1)
          ver(2)=z(2)*x(1)-x(2)*z(1)
          ver(3)=x(2)*y(1)-y(2)*x(1)
          DO i=1,3
              vrr(i)=0.0D0
              DO j=1,3
                  vrr(i)=vrr(i)+d0(i,j)*ver(j)
              END DO
          END DO
          norm=0.0D0
          DO i=1,3
              norm=norm+vrr(i)**2
          END DO
          norm=DSQRT(norm)
          DO i=1,3
              vrr(i)=vrr(i)/norm
          END DO
          alpha=DACOS(vrr(1)*pax(1)+vrr(2)*pax(2)+vrr(3)*pax(3))

*=======================================================================
*---- Compute the sign of the angle ALPHA ------------------------------
*=======================================================================

          vsg(1)=vrr(2)*pax(3)-pax(2)*vrr(3)
          vsg(2)=vrr(3)*pax(1)-pax(3)*vrr(1)
          vsg(3)=vrr(1)*pax(2)-pax(1)*vrr(2)
          norm=0.0D0
          DO i=1,3
              norm=norm+vsg(i)**2
          END DO       
          norm=DSQRT(norm)
          DO i=1,3
              vsg(i)=vsg(i)/norm
          END DO
          norm=pax(1)*ax(1)+pax(2)*ax(2)+pax(3)*ax(3)
          sgns=vsg(1)*ax(1)+vsg(2)*ax(2)+vsg(3)*ax(3)
          sgns=DSIGN(1.0D0,sgns)
          alpha=(0.50D0*pi-alpha)*sgns
          CALL comd0(a,b,c,alpha,d0)
      END IF

*=======================================================================
*---- Finally rotate the primitive coordinates and find the ------------
*---- required hydrogen coordinates ------------------------------------
*=======================================================================

      DO m=1,2
          x0(m)=d0(1,1)*x(m)+d0(1,2)*y(m)+d0(1,3)*z(m)
          y0(m)=d0(2,1)*x(m)+d0(2,2)*y(m)+d0(2,3)*z(m)
          z0(m)=d0(3,1)*x(m)+d0(3,2)*y(m)+d0(3,3)*z(m)
      END DO
      DO m=1,2
          x0(m)=x0(m)*length
          y0(m)=y0(m)*length
          z0(m)=z0(m)*length
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
