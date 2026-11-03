      subroutine meta_bond(at1,at2,ntap,xp0,yp0,zp0,co,oc,rc,ddx,ddy
     &     ,ddz)
      implicit none
      integer i
      integer ntap
      integer at1, at2
      real*8 rc
      real*8 dx0,dy0,dz0,dx1,dy1,dz1
      real*8 co(3,3),oc(3,3)
      real*8 xp0(*),yp0(*),zp0(*)
      real*8 xp1(ntap),yp1(ntap),zp1(ntap)
      real*8 ddx(*),ddy(*),ddz(*)

      include 'pbc.h'
      
      do i = 1,ntap
        ddx(i) = 0.d0
        ddy(i) = 0.d0
        ddz(i) = 0.d0
      enddo
      
c change frame
      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xp1,yp1,zp1) 
      
      dx1 = xp1(at2) - xp1(at1) 
      dy1 = yp1(at2) - yp1(at1) 
      dz1 = zp1(at2) - zp1(at1) 
      
c apply pbc
      
      dx1 = dx1 - 2.d0 * pbc(dx1)
      dy1 = dy1 - 2.d0 * pbc(dy1)
      dz1 = dz1 - 2.d0 * pbc(dz1)
      
c back to cartesian
      
      dx0=co(1,1)*dx1+co(1,2)*dy1+co(1,3)*dz1
      dy0=co(2,1)*dx1+co(2,2)*dy1+co(2,3)*dz1
      dz0=co(3,1)*dx1+co(3,2)*dy1+co(3,3)*dz1
      
      rc = dsqrt(dx0*dx0+dy0*dy0+dz0*dz0)
      
      ddx(at2) = dx0/rc
      ddy(at2) = dy0/rc
      ddz(at2) = dz0/rc
      ddx(at1) = -dx0/rc
      ddy(at1) = -dy0/rc
      ddz(at1) = -dz0/rc
      
      
      return
      end

      SUBROUTINE meta_bend(at1,at2,at3,ntap,xp0,yp0,zp0,co,oc,rc,ddx,ddy
     &     ,ddz)
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------
      
      INTEGER at1,at2,at3,at4,ntap
      REAL*8  rc,xp0(*),yp0(*),zp0(*),co(3,3),oc(3,3),ddx(*),ddy(*),ddz(
     &     *)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i
      REAL*8  xr1,xr2,xr3,yr1,yr2,yr3,zr1,zr2,zr3,x12,x32,y12,y32,z12
     &     ,z32,rs12,rs32,uux1,uux2,uux3,uuy1,uuy2,uuy3,uuz1,uuz2,uuz3
     &     ,xr31,yr31,zr31,rs31,rsp31,x31,y31,z31,rsp12,rsp32,k12,r1,r2
      REAL*8 x12c,y12c,z12c,x32c,y32c,z32c
      REAL*8  dcc2,cb,sb,bb,pi,ssb
      real*8 xp1(ntap),yp1(ntap),zp1(ntap)
!      include 'pbc.h'

*==================== EXECUTABLE STATEMENTS ============================

      do i = 1,ntap
        ddx(i) = 0.d0
        ddy(i) = 0.d0
        ddz(i) = 0.d0
      enddo

      pi=4.0D0*DATAN(1.0D0)

      xr1=xp0(at1)
      yr1=yp0(at1)
      zr1=zp0(at1)
      xr2=xp0(at2)
      yr2=yp0(at2)
      zr2=zp0(at2)
      xr3=xp0(at3)
      yr3=yp0(at3)
      zr3=zp0(at3)
      x12=xr1-xr2
      y12=yr1-yr2
      z12=zr1-zr2
      x32=xr3-xr2
      y32=yr3-yr2
      z32=zr3-zr2

c change frame
!      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xp1,yp1,zp1) 

!      x12c = xp1(at1) - xp1(at2)
!      y12c = yp1(at1) - yp1(at2)
!      z12c = zp1(at1) - zp1(at2)
!      x32c = xp1(at3) - xp1(at2)
!      y32c = yp1(at3) - yp1(at2)
!      z32c = zp1(at3) - zp1(at2)

c apply pbc
      
!      x12c = x12c - 2.d0 * pbc(x12c)
!      y12c = y12c - 2.d0 * pbc(y12c)
!      z12c = z12c - 2.d0 * pbc(z12c)
!      x32c = x32c - 2.d0 * pbc(x32c)
!      y32c = y32c - 2.d0 * pbc(y32c)
!      z32c = z32c - 2.d0 * pbc(z32c)
      
c back to cartesian
      
!      x12=co(1,1)*x12c+co(1,2)*y12c+co(1,3)*z12c
!      y12=co(2,1)*x12c+co(2,2)*y12c+co(2,3)*z12c
!      z12=co(3,1)*x12c+co(3,2)*y12c+co(3,3)*z12c
!      x32=co(1,1)*x32c+co(1,2)*y32c+co(1,3)*z32c
!      y32=co(2,1)*x32c+co(2,2)*y32c+co(2,3)*z32c
!      z32=co(3,1)*x32c+co(3,2)*y32c+co(3,3)*z32c

      rs12=x12**2+y12**2+z12**2
      rs32=x32**2+y32**2+z32**2
      
      dcc2=DSQRT(rs12*rs32)
      cb=(x12*x32+y12*y32+z12*z32)/dcc2
      sb=DSQRT(1.0d0-cb**2)
      bb=DACOS(cb)
      rc = bb*180.0d0/pi
      
      ssb = dsqrt((y12*z32 - y12*z12)**2 + (z12*x32 - z32*x12)**2 + 
     &     (x12*y32 - x32*y12)**2) / dcc2
      
      if(ssb.lt.0.0) rc = - rc
      
      uux1=x32/dcc2-cb*x12/rs12
      uux3=x12/dcc2-cb*x32/rs32
      uux2=-uux1-uux3
      uuy1=y32/dcc2-cb*y12/rs12
      uuy3=y12/dcc2-cb*y32/rs32
      uuy2=-uuy1-uuy3
      uuz1=z32/dcc2-cb*z12/rs12
      uuz3=z12/dcc2-cb*z32/rs32
      uuz2=-uuz1-uuz3
      ddx(at1)=-uux1*180.0/pi
      ddx(at2)=-uux2*180.0/pi
      ddx(at3)=-uux3*180.0/pi
      ddy(at1)=-uuy1*180.0/pi
      ddy(at2)=-uuy2*180.0/pi
      ddy(at3)=-uuy3*180.0/pi
      ddz(at1)=-uuz1*180.0/pi
      ddz(at2)=-uuz2*180.0/pi
      ddz(at3)=-uuz3*180.0/pi
      
*================= END OF EXECUTABLE STATEMENTS ========================
      
      RETURN
      END
      
      
      SUBROUTINE meta_tors(at1,at2,at3,at4,ntap,xp0,yp0,zp0,co,oc,rc,ddx
     &     ,ddy,ddz)
      
      IMPLICIT none
      
*----------------------- ARGUMENTS -------------------------------------
      
      INTEGER at1,at2,at3,at4,ntap
      REAL*8  rc,xp0(*),yp0(*),zp0(*),co(3,3),oc(3,3),ddx(*),ddy(*),ddz(
     &     *)
      
*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j
      real*8  x12,y12,z12,r12,dfunc
      REAL*8  xr1,yr1,zr1,xr2,yr2,zr2,pi
      REAL*8  xr3,yr3,zr3,xr4,yr4,zr4,x21,y21,z21,x32,y32,z32,x43,y43
     &     ,z43
      REAL*8  nx1,ny1,nz1,nx2,ny2,nz2,n1mod,n2mod,aaa,bbb,ccc,ddd,sign
      REAL*8  cosdrd,sindrd,angle,drd,uux1,uuy1,uuz1,uux2,uuy2,uuz2
      REAL*8  uux3,uuy3,uuz3,uux4,uuy4,uuz4,rsp32,rsq32

      do i = 1,ntap
        ddx(i) = 0.d0
        ddy(i) = 0.d0
        ddz(i) = 0.d0
      enddo

      pi=dacos(-1.0d0)      
      
      xr1=xp0(at1)
      yr1=yp0(at1)
      zr1=zp0(at1)
      xr2=xp0(at2)
      yr2=yp0(at2)
      zr2=zp0(at2)
      xr3=xp0(at3)
      yr3=yp0(at3)
      zr3=zp0(at3)
      xr4=xp0(at4)
      yr4=yp0(at4)
      zr4=zp0(at4)
      x21=xr2-xr1
      y21=yr2-yr1
      z21=zr2-zr1
      x32=xr3-xr2
      y32=yr3-yr2
      z32=zr3-zr2
      x43=xr4-xr3
      y43=yr4-yr3
      z43=zr4-zr3
      rsq32=x32**2+y32**2+z32**2
      rsp32=DSQRT(rsq32)

      nx1=y21*z32-y32*z21
      ny1=z21*x32-z32*x21
      nz1=x21*y32-x32*y21

      nx2=y32*z43-y43*z32
      ny2=z32*x43-z43*x32
      nz2=x32*y43-x43*y32

      n1mod=(nx1*nx1)+(ny1*ny1)+(nz1*nz1)
      n1mod=sqrt(n1mod)
      n2mod=(nx2*nx2)+(ny2*ny2)+(nz2*nz2)
      n2mod=sqrt(n2mod)
      aaa=(nx1*nx2)+(ny1*ny2)+(nz1*nz2)
      bbb=n1mod*n2mod
      cosdrd=aaa/bbb
      
      ccc=(x32*(ny2*nz1-nz2*ny1)+y32*(nx1*nz2-nz1*nx2)+
     &     z32*(nx2*ny1-ny2*nx1))
      ddd=n2mod*n1mod*rsp32
      sindrd=ccc/ddd
      drd=atan2(sindrd,cosdrd)
      angle = -(drd*(180.))/(pi)

      if(angle.gt.180.) angle=180.0d0
      if(angle.lt.-180.) angle=-180.0d0
      
      rc = angle

      if(sindrd.ge.0.0d0) then
         sign=+1.0d0
      else
         sign=-1.0d0
      end if
      
      dfunc=sign*(-180.d0/(pi))*(-1.0d0/dsqrt(1.0d0-cosdrd**2))
      
      uux1 = (z32*ny2-y32*nz2) -aaa*(z32*ny1-y32*nz1)/n1mod**2
      uux2 = (-(z21+z32)*ny2 + z43*ny1 + (y32+y21)*nz2 - y43*nz1)
     &     -aaa*((-ny1*(z21+z32)+nz1*(y32+y21))/n1mod**2+(z43*ny2-y43
     &     *nz2)/n2mod**2)
      uux3 = (z21*ny2 - (z32+z43)*ny1 - y21*nz2 + (y43+y32)*nz1)
     &     -aaa*((z21*ny1-y21*nz1)/n1mod**2+((y43+y32)*nz2-(z32+z43)*ny2
     &     )/n2mod**2)
      uux4 = (z32*ny1-y32*nz1) - aaa*(z32*ny2-y32*nz2)/n2mod**2
      
      
      
      uuy1 = (-z32*nx2+x32*nz2) -aaa*(-z32*nx1+x32*nz1)/n1mod**2
      uuy2 = ((z21+z32)*nx2 - z43*nx1 - (x32+x21)*nz2 + x43*nz1)
     &     -aaa*((nx1*(z21+z32)-nz1*(x32+x21))/n1mod**2+(-z43*nx2+x43
     &     *nz2)/n2mod**2)
      uuy3 = (-z21*nx2 + (z32+z43)*nx1 + x21*nz2 - (x43+x32)*nz1)
     &     -aaa*((-z21*nx1+x21*nz1)/n1mod**2+(-(x43+x32)*nz2+(z32+z43)
     &     *nx2)/n2mod**2)
      uuy4 = (-z32*nx1+x32*nz1) - aaa*(-z32*nx2+x32*nz2)/n2mod**2
      
      
      
      uuz1 = (-x32*ny2+y32*nx2) -aaa*(-x32*ny1+y32*nx1)/n1mod**2
      uuz2 = ((x21+x32)*ny2 - x43*ny1 - (y32+y21)*nx2 + y43*nx1)
     &     -aaa*((ny1*(x21+x32)-nx1*(y32+y21))/n1mod**2+(-x43*ny2+y43
     &     *nx2)/n2mod**2)
      uuz3 = (-x21*ny2 + (x32+x43)*ny1 + y21*nx2 - (y43+y32)*nx1)
     &     -aaa*((-x21*ny1+y21*nx1)/n1mod**2+(-(y43+y32)*nx2+(x32+x43)
     &     *ny2)/n2mod**2)
      uuz4 = (-x32*ny1+y32*nx1) - aaa*(-x32*ny2+y32*nx2)/n2mod**2
      
      ddx(at1)=dfunc*uux1/bbb
      ddx(at2)=dfunc*uux2/bbb
      ddx(at3)=dfunc*uux3/bbb
      ddx(at4)=dfunc*uux4/bbb
      ddy(at1)=dfunc*uuy1/bbb
      ddy(at2)=dfunc*uuy2/bbb
      ddy(at3)=dfunc*uuy3/bbb
      ddy(at4)=dfunc*uuy4/bbb
      ddz(at1)=dfunc*uuz1/bbb
      ddz(at2)=dfunc*uuz2/bbb
      ddz(at3)=dfunc*uuz3/bbb
      ddz(at4)=dfunc*uuz4/bbb
      
      return
      end
      
