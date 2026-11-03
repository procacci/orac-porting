c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
       real*8 function dist_ijk(ni,nj,nk,dx,dy,dz,co)

       implicit none
c----------------------------
c declarations
c----------------------------
      real*8 dx,dy,dz,co(3,3)
      integer ni,nj,nk

      real*8 d,dmin,dt
      real*8 lx,ly,lz
      real*8 mx,my,mz
      real*8 dmx,dmy,dmz

      integer nv(10,3)
      integer i,j,imin,jmin
      integer ndtmax

c----------------------------
c initialisations
c----------------------------
      ndtmax=500

      nv(1,1)=0
      nv(2,1)=0
      nv(3,1)=0
      nv(4,1)=0
      nv(5,1)=1
      nv(6,1)=1
      nv(7,1)=1
      nv(8,1)=1

      nv(1,2)=0
      nv(2,2)=0
      nv(3,2)=1
      nv(4,2)=1
      nv(5,2)=0
      nv(6,2)=0
      nv(7,2)=1
      nv(8,2)=1

      nv(1,3)=0
      nv(2,3)=1
      nv(3,3)=0
      nv(4,3)=1
      nv(5,3)=0
      nv(6,3)=1
      nv(7,3)=0
      nv(8,3)=1

c----------------------------
c calcul de la distance
c minimale entre les coins
c des cellules (0,0,0)
c et (ni,nj,nk)
c----------------------------
      dmin=1.d+08

      do i=1,8
         do j=1,8

            lx=(ni+nv(j,1)-nv(i,1))*dx
            ly=(nj+nv(j,2)-nv(i,2))*dy
            lz=(nk+nv(j,3)-nv(i,3))*dz

            mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
            my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
            mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

            d=mx*mx+my*my+mz*mz

            if(d.lt.dmin) then
               dmin=d
               imin=i
               jmin=j
            endif

         end do
      end do

c----------------------------
c on fait varier la distance
c trouvee precedement pour
c voir si la distance 
c minimale ne se trouve pas
c sur une arrete
c----------------------------
      lx=(ni+nv(jmin,1)-nv(imin,1))*dx
      ly=(nj+nv(jmin,2)-nv(imin,2))*dy
      lz=(nk+nv(jmin,3)-nv(imin,3))*dz

      mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
      my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
      mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

c----------------------------
c sur l'axe x
c----------------------------
      dt=sign(1.,0.5-nv(jmin,1))*dx/ndtmax

      dmin=1.d+08

      do i=1,ndtmax
         dmx=i*co(1,1)*dt
         dmy=i*co(2,1)*dt
         dmz=i*co(3,1)*dt
         d=(mx+dmx)**2+(my+dmy)**2+(mz+dmz)**2
         if(d.lt.dmin) dmin=d
      end do

c----------------------------
c sur l'axe y
c----------------------------
      dt=sign(1.,0.5-nv(jmin,2))*dy/ndtmax

      do i=1,ndtmax
         dmx=i*co(1,2)*dt
         dmy=i*co(2,2)*dt
         dmz=i*co(3,2)*dt
         d=(mx+dmx)**2+(my+dmy)**2+(mz+dmz)**2
         if(d.lt.dmin) dmin=d
      end do

c----------------------------
c sur l'axe z
c----------------------------
      dt=sign(1.,0.5-nv(jmin,3))*dz/ndtmax

      do i=1,ndtmax
         dmx=i*co(1,3)*dt
         dmy=i*co(2,3)*dt
         dmz=i*co(3,3)*dt
         d=(mx+dmx)**2+(my+dmy)**2+(mz+dmz)**2
         if(d.lt.dmin) dmin=d
      end do

      dist_ijk=dmin
      return
      end
