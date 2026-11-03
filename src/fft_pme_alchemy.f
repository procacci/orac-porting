      subroutine fft_pme_alchemy(ntap,xp,yp,zp,xpcm,ypcm,zpcm,pmechg
     &     ,charg0,co,oc,volume,alphal,order,nfft1,nfft2,nfft3,eer,fpx
     &     ,fpy,fpz,stressc,atomp,grppt,pressure,rkcut,dwrk)
      
      use spme
      
      IMPLICIT none

c---  stuff for b-spline interpolation and FFT 

      INTEGER   NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER,i,j

      REAL*8  theta1(mth),dtheta1(mth),theta2(mth),dtheta2(mth)
     &     ,theta3(mth),dtheta3(mth),rkcut
      REAL*8  Q(maxt),fr1(max_atm),fr2(max_atm),fr3(max_atm),recip(3,3)
     &     ,x(max_atm),y(max_atm),z(max_atm),fx(max_atm),fy(max_atm)
     &     ,fz(max_atm),sgsx,sgsy,sgsz
     &     ,vir(3,3)
      REAL*8 sx,sy,sz,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9
      INTEGER k

c---  non local stuff 
      
      INTEGER  ntap,atomp(*),grppt(2,*),count,mia,ia
      REAL*8   fpx(*),fpy(*),fpz(*)
      REAL*8   xp(*),yp(*),zp(*),xpcm(*),ypcm(*),zpcm(*),pmechg(*)
     &     ,charg0(*),dwrk
      REAL*8   co(3,3),oc(3,3),volume,alphal,stressc(3,3),eer
      LOGICAL pressure

c---  transfor all to cartesian coordinates

      
      numatoms=ntap
      do i=1,numatoms
         x(i)=co(1,1)*xp(i)+co(1,2)*yp(i)+co(1,3)*zp(i)
         y(i)=co(2,1)*xp(i)+co(2,2)*yp(i)+co(2,3)*zp(i)
         z(i)=co(3,1)*xp(i)+co(3,2)*yp(i)+co(3,3)*zp(i)
      end do
      do i=1,numatoms
         fx(i)=0.0D0
         fy(i)=0.0D0
         fz(i)=0.0D0
      end do

c--   call stand-alone Darden routine 

      do i=1,3
         do j=1,3
            recip(i,j) = 0.5*oc(j,i)
         end do
      end do

#ifdef FFTW
      call do_pmesh_alchemy(numatoms,x,y,z,pmechg,charg0,recip,volume
     &     ,alphal,order,nfft1,nfft2,nfft3,eer,fx,fy,fz,vir,bsp_mod1
     &     ,bsp_mod2,bsp_mod3,theta1,theta2,theta3,dtheta1,dtheta2
     &     ,dtheta3,fr1,fr2,fr3,rkcut,planf,planr,dwrk)
#else
      call do_pmesh_alchemy(numatoms,x,y,z,pmechg,charg0,recip,volume
     &     ,alphal,order,nfft1,nfft2,nfft3,eer,fx,fy,fz,vir,sizfftab
     &     ,sizffwrk,siztheta,siz_Q,bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q
     &     ,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fr1,fr2
     &     ,fr3,rkcut,dwrk)
#endif

      IF(pressure) THEN
         stc1=0.0D0
         stc2=0.0D0
         stc3=0.0D0
         stc4=0.0D0
         stc5=0.0D0
         stc6=0.0D0
         stc7=0.0D0
         stc8=0.0D0
         stc9=0.0D0
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE)
!$OMP& SHARED(xpcm,ypcm,zpcm,co,fx,fy,fz,xp,yp,zp,atomp,ntap)
!$OMP& REDUCTION(+:stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9)
         DO i=1,ntap
            j=atomp(i)
            sx=xp(i)-xpcm(j)
            sy=yp(i)-ypcm(j)
            sz=zp(i)-zpcm(j)
            sgsx=sx*co(1,1)+sy*co(1,2)+sz*co(1,3)
            sgsy=sx*co(2,1)+sy*co(2,2)+sz*co(2,3)
            sgsz=sx*co(3,1)+sy*co(3,2)+sz*co(3,3)
            stc1=stc1-sgsx*fx(i)
            stc2=stc2-sgsy*fx(i)
            stc3=stc3-sgsz*fx(i)
            stc4=stc4-sgsx*fy(i)
            stc5=stc5-sgsy*fy(i)
            stc6=stc6-sgsz*fy(i)
            stc7=stc7-sgsx*fz(i)
            stc8=stc8-sgsy*fz(i)
            stc9=stc9-sgsz*fz(i)
         END DO
!$OMP END PARALLEL DO 
         stressc(1,1)=stressc(1,1)-vir(1,1)+stc1
         stressc(1,2)=stressc(1,2)-vir(1,2)+stc2
         stressc(1,3)=stressc(1,3)-vir(1,3)+stc3
         stressc(2,1)=stressc(2,1)-vir(2,1)+stc4
         stressc(2,2)=stressc(2,2)-vir(2,2)+stc5
         stressc(2,3)=stressc(2,3)-vir(2,3)+stc6
         stressc(3,1)=stressc(3,1)-vir(3,1)+stc7
         stressc(3,2)=stressc(3,2)-vir(3,2)+stc8
         stressc(3,3)=stressc(3,3)-vir(3,3)+stc9
      END IF

      do i=1,numatoms
         fpx(i) = fpx(i)+fx(i)
         fpy(i) = fpy(i)+fy(i)
         fpz(i) = fpz(i)+fz(i)
      end do

      return
      end 
