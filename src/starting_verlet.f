      SUBROUTINE starting_verlet(mass,ntap,xp0,yp0,zp0,vpx,vpy,vpz,fpxn
     &     ,fpyn,fpzn,dt,max_dist)

************************************************************************
*                                                                      *
*  STRARTING_VERLET advance coordinates for one time step and          *
*  velocities for half time step. This routine is used for             *
*  poorly equilibrated samples with chane for bad contacts             *
*  Coordinates can never be advanced more than a prefixed value given  *
*  as MAX_DIST                                                         *
*  Coordinates and velocities arrays are overwritten.                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      integer      ntap
      real*8       dt,mass(*)
      REAL*8       xp0(*),yp0(*),zp0(*),max_dist(3)
      REAL*8       vpx(*),vpy(*),vpz(*)
      REAL*8       fpxn(*),fpyn(*),fpzn(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8       ts2,tfact,dx,dy,dz
      INTEGER      i

*==================== EXECUTABLE STATEMENTS ============================

      ts2=dt*dt

*=======================================================================
*  advances atomic coordinates for a full time step 
*  and advances solute velocities for half time step 
*=======================================================================

      do i=1,ntap
         tfact=0.5*ts2/mass(i)
         dx = dt*vpx(i)+tfact*fpxn(i)
         dy = dt*vpy(i)+tfact*fpyn(i)
         dz = dt*vpz(i)+tfact*fpzn(i)
         if(abs(dx).gt.max_dist(1)) dx = sign(1.d0,dx)*max_dist(1)  
         if(abs(dy).gt.max_dist(2)) dy = sign(1.d0,dy)*max_dist(2)  
         if(abs(dz).gt.max_dist(3)) dz = sign(1.d0,dz)*max_dist(3)  
         xp0(i)=xp0(i)+dx
         yp0(i)=yp0(i)+dy
         zp0(i)=zp0(i)+dz
         tfact=0.5*dt/mass(i)
         vpx(i)=vpx(i) + tfact*fpxn(i)
         vpy(i)=vpy(i) + tfact*fpyn(i)
         vpz(i)=vpz(i) + tfact*fpzn(i)
      end do   

      RETURN
      END

