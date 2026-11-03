#ifdef FFTW 
subroutine fft_pme_init(numatoms,nfft1,nfft2,nfft3,order,bsp_mod1,bsp_mod2,bsp_mod3,planf,planb)
  use, intrinsic :: iso_c_binding
#ifdef _OMP_
  use omp_integr, only: nthr2
#endif
  IMPLICIT NONE
  include 'fftw3.f03'  
#else
  subroutine fft_pme_init(numatoms,nfft1,nfft2,nfft3,order,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork)
    IMPLICIT NONE
#endif
!***********************************************************************
!
!     To be called from MTSMD
!---  ON INPUT
!     nfft1,nfft2,nfft3: grid points in the k1,k2,k3 directions
!     order            : order of B-spline interpolation
!---  ON OUTPUT 
!     sizfftab is permanent 3d fft table storage
!     sizffwrk is temporary 3d fft work storage
!     siztheta is size of arrays theta1-3 dtheta1-3
!     sizheap is total size of permanent storage
!     sizstack is total size of temporary storage
!     bsp_mod1-3 hold the moduli of the inverse DFT of the B splines
!
!***********************************************************************

INTEGER    NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER
#ifndef FFTW
INTEGER   sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
REAL*8 bsp_mod1(*),bsp_mod2(*),bsp_mod3(*)
REAL*8 fftable(*),ffwork(*)
#endif

#ifdef FFTW
!     initialization of work array is no longer needed if fftw is used. 
REAL*8 bsp_mod1(*),bsp_mod2(*),bsp_mod3(*)
real*8, allocatable ::  Q(:,:,:)
complex*16, allocatable :: FQ(:,:,:)
integer ierr
type(C_PTR) :: planf,planb

#ifdef _OMP_
integer void
void=fftw_init_threads()
call  fftw_plan_with_nthreads(nthr2)
#endif

allocate(Q(nfft1,nfft2,nfft3),stat=ierr) 
allocate(FQ(1+nfft1/2,nfft2,nfft3),stat=ierr) 

planf=fftw_plan_dft_r2c_3d(nfft3,nfft2,nfft1,Q,FQ,FFTW_ESTIMATE)
planb=fftw_plan_dft_c2r_3d(nfft3,nfft2,nfft1,FQ,Q,FFTW_ESTIMATE)
call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,order)
deallocate(Q,FQ,stat=ierr)
#else
call pmesh_kspace_get_sizes(nfft1,nfft2,nfft3,numatoms,order,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
call pmesh_kspace_setup(bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
#endif
return
end subroutine fft_pme_init

#ifdef FFTW
subroutine set_vir_scalar_sum(recip,nfft1,nfft2,nfft3,face,fac12,fac13,fac23)
  implicit none
  integer nfft1,nfft2,nfft3
  REAL*8 recip(3,3),face(*),mhat1,mhat2,mhat3
  integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
  integer nf1,nf2,nf3,ndim1
  logical fac12(*),fac13(*),fac23(*)
  ndim1=1+nfft1/2
  indtop = ndim1*nfft2*nfft3
  nff = ndim1*nfft2
  nf1 = nfft1/2
  if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
  nf2 = nfft2/2
  if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
  nf3 = nfft3/2
  if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
  do ind = 1,indtop-1
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
     if(k1.eq.1.or.(k1.eq.ndim1.and.mod(nfft1,2).eq.0)) THEN 
        face(ind)=1.d0
     else
        face(ind)=2.0
     endif

     if(mhat1*mhat2.lt.0.d0.and.k2.le.ndim1) THEN
        fac12(ind)=.false.
     else
        fac12(ind)=.true.
     end if
     if(mhat1*mhat3.lt.0.d0.and.k3.le.ndim1) THEN
        fac13(ind)=.false.
     else
        fac13(ind)=.true.
     end if
     
     if(k2.le.ndim1.and.k3.le.ndim1) THEN
        if(mhat2*mhat3.lt.0.d0)  THEN
           fac23(ind)=.false.
        else
           fac23(ind)=.true.
        endif
     else if(k2.gt.ndim1.and.k3.le.ndim1) THEN
        if(mhat2*mhat3.lt.0.d0)  THEN
           fac23(ind)=.true.
        else
           fac23(ind)=.false.
        endif
     else if(k3.gt.ndim1.and.k2.le.ndim1) THEN
        if(mhat2*mhat3.lt.0.d0)  THEN
           fac23(ind)=.true.
        else
           fac23(ind)=.false.
        endif
     else if(k3.gt.ndim1.and.k2.gt.ndim1) THEN
        if(mhat2*mhat3.lt.0.d0)  THEN
           fac23(ind)=.false.
        else
           fac23(ind)=.true.
        endif
     end if
  end do
  return
end subroutine set_vir_scalar_sum
#endif
