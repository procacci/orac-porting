#ifdef FFTW
      subroutine scalar_sum(
     $     Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1
     &     ,nfft2,nfft3,eer,vir,rkcut)
      implicit none
      integer nfft1,nfft2,nfft3
      complex*16 Q(1+nfft1/2,nfft2,nfft3)
      REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume,rkcut
      REAL*8 eer,vir(3,3)
      REAL*8 recip(3,3)
      REAL*8 pi,fac,denom,eterm,vterm,energy
      integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3,ndim1
      REAL*8 mhat1,mhat2,mhat3,msq,struc2,rkcut2,fac12,fac13,fac23,vvir
      ndim1=1+nfft1/2
      indtop = ndim1*nfft2*nfft3
      pi = 3.14159265358979323846
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/ewaldcof**2
      nff = ndim1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      energy = 0.d0
      DO k1 = 1,3
         DO k2 = 1,3
            vir(k1,k2) = 0.0D0
         END DO
      END DO
#ifdef _OMP_
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC)  
!$OMP& SHARED(nfft1,nfft2,nfft3,recip,bsp_mod1,bsp_mod2,bsp_mod3)  
!$OMP& SHARED(pi,volume,Q,ndim1,nff,nf1,nf2,nf3,rkcut2,fac,indtop) 
!$OMP& REDUCTION(+:energy,vir)  
#endif
      do 100 ind = 1,indtop-1
c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1
         k3 = ind/nff + 1
         jnd = ind - (k3-1)*nff
         k2 = jnd/ndim1 + 1
         k1 = jnd - (k2-1)*ndim1 +1
         m1 = k1 - 1
         if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
         m2 = k2 - 1
         if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
         m3 = k3 - 1
         if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
         mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
         mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
         mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
         msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
         IF(msq .GT. rkcut2) THEN 
            eterm = 0.d0 
         ELSE
            denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
            eterm = dexp(-fac*msq)/denom
            vterm = 2.d0*(fac*msq + 1.d0)/msq
            if(k1.eq.1.or.(k1.eq.ndim1.and.mod(nfft1,2).eq.0)) THEN 
              struc2 =      eterm*Q(k1,k2,k3)*dconjg(Q(k1,k2,k3))
              vvir=struc2*vterm
              fac12=1.d0
              fac13=1.d0
              fac23=1.d0
            ELSE
              struc2 = 2.d0*eterm*Q(k1,k2,k3)*dconjg(Q(k1,k2,k3))
              vvir=struc2*vterm
              if(mhat1*mhat2.lt.0.d0.and.k2.le.ndim1) THEN
                fac12=0.d0
              else
                fac12=1.d0
              end if
              if(mhat1*mhat3.lt.0.d0.and.k3.le.ndim1) THEN
                fac13=0.d0
              else
                fac13=1.d0
              end if

              if(k2.le.ndim1.and.k3.le.ndim1) THEN
                if(mhat2*mhat3.lt.0.d0)  THEN
                  fac23=0.d0
                else
                  fac23=1.d0
                endif
              else if(k2.gt.ndim1.and.k3.le.ndim1) THEN
                if(mhat2*mhat3.lt.0.d0)  THEN
                  fac23=1.d0
                else
                  fac23=0.d0
                endif
              else if(k3.gt.ndim1.and.k2.le.ndim1) THEN
                if(mhat2*mhat3.lt.0.d0)  THEN
                  fac23=1.d0
                else
                  fac23=0.d0
                endif
              else if(k3.gt.ndim1.and.k2.gt.ndim1) THEN
                if(mhat2*mhat3.lt.0.d0)  THEN
                  fac23=0.d0
                else
                  fac23=1.d0
                endif
              end if
            END IF
            energy = energy + struc2
            vir(1,1) = vir(1,1) + struc2*(vterm*mhat1*mhat1 -1.d0)
            vir(2,2) = vir(2,2) + struc2*(vterm*mhat2*mhat2 -1.d0)
            vir(3,3) = vir(3,3) + struc2*(vterm*mhat3*mhat3 -1.d0)
            vir(1,2) = vir(1,2) + vvir*mhat1*mhat2*fac12
            vir(1,3) = vir(1,3) + vvir*mhat1*mhat3*fac13
            vir(2,3) = vir(2,3) + vvir*mhat2*mhat3*fac23
         END IF
!        B,C and F[Q] are in reciprocal space; H_rec and Q are in direct space.
!        this BC*F[Q] = F[H_rec]*F[Q] where H_rec=F^-1[BC]  
!        so that F^-1(F[H]*F(Q)) =H_rec x Q (where "x" indicate convolution) 
1098     Q(k1,k2,k3) = eterm * Q(k1,k2,k3) ! array is overwritten by the product BC*F[Q]
100   continue
#ifdef _OMP_
!$OMP END PARALLEL DO
#endif
      eer = 0.5d0 * energy
      vir(2,1)=vir(1,2)
      vir(3,1)=vir(1,3)
      vir(3,2)=vir(2,3)
      DO k1 = 1,3
         DO k2 = 1,3
            vir(k1,k2) = 0.5d0*vir(k1,k2)
         END DO
      END DO
      return
      end
#else
      subroutine scalar_sum(
     $         Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,vir,rkcut)
      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3)
      REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume,rkcut
      REAL*8 eer,vir(3,3)
      REAL*8 recip(3,3)
      REAL*8 pi,fac,denom,eterm,vterm,energy
      integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3
      REAL*8 mhat1,mhat2,mhat3,msq,struc2,rkcut2
      indtop = nfft1*nfft2*nfft3
      pi = 3.14159265358979323846
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/ewaldcof**2
      nff = nfft1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      energy = 0.d0
      DO k1 = 1,3
         DO k2 = 1,3
            vir(k1,k2) = 0.0D0
         END DO
      END DO
      do 100 ind = 1,indtop-1
c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1
         k3 = ind/nff + 1
         jnd = ind - (k3-1)*nff
         k2 = jnd/nfft1 + 1
         k1 = jnd - (k2-1)*nfft1 +1
         m1 = k1 - 1
         if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
         m2 = k2 - 1
         if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
         m3 = k3 - 1
         if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
         mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
         mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
         mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
         msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
         IF(msq .GT. rkcut2) THEN 
            eterm = 0.d0 
         ELSE
            denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
            eterm = dexp(-fac*msq)/denom
            vterm = 2.d0*(fac*msq + 1.d0)/msq
            struc2 = Q(1,k1,k2,k3)**2 + Q(2,k1,k2,k3)**2
            energy = energy + eterm * struc2
            vir(1,1) = vir(1,1) + eterm * struc2 * (vterm*mhat1*mhat1 -
     &           1.d0)
            vir(1,2) = vir(1,2) + eterm * struc2 * (vterm*mhat1*mhat2)
            vir(1,3) = vir(1,3) + eterm * struc2 * (vterm*mhat1*mhat3)
            vir(2,2) = vir(2,2) + eterm * struc2 * (vterm*mhat2*mhat2 -
     &           1.d0)
            vir(2,3) = vir(2,3) + eterm * struc2 * (vterm*mhat2*mhat3)
            vir(3,3) = vir(3,3) + eterm * struc2 * (vterm*mhat3*mhat3 -
     &           1.d0)
         END IF
1098     Q(1,k1,k2,k3) = eterm * Q(1,k1,k2,k3)
         Q(2,k1,k2,k3) = eterm * Q(2,k1,k2,k3)
100   continue
      eer = 0.5d0 * energy
      vir(2,1)=vir(1,2)
      vir(3,1)=vir(1,3)
      vir(3,2)=vir(2,3)
      DO k1 = 1,3
         DO k2 = 1,3
            vir(k1,k2) = 0.5d0*vir(k1,k2)
         END DO
      END DO
      return
      end
#endif
