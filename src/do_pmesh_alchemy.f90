#ifdef FFTW
subroutine do_pmesh_alchemy(numatoms,x,y,z,charge,charge0,recip,volume,ewald_coeff,order, & 
     nfft1,nfft2,nfft3,eer,dx,dy,dz,virial,bsp_mod1,bsp_mod2,bsp_mod3,theta1,theta2,theta3, &
     dtheta1,dtheta2,dtheta3,fr1,fr2,fr3,rkcut,planf,planr,dwrk)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'  
#else
subroutine do_pmesh_alchemy(numatoms,x,y,z,charge,charge0,recip,volume,ewald_coeff, & 
     order,nfft1,nfft2,nfft3, eer,dx,dy,dz,virial,sizfftab,sizffwrk,siztheta,siz_Q, & 
     bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3, & 
     fr1,fr2,fr3,rkcut,dwrk)
  implicit none
#endif
!-------------------------------------------------------------------------------------------
! INPUT 
  !       numatoms:  number of atoms
!       x,y,z:   atomic coords
!       charge  atomic charges
!       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
!       volume: the volume of the unit cell
!       ewald_coeff:   ewald convergence parameter
!       order: the order of Bspline interpolation. E.g. cubic is order 4
!          fifth degree is order 6 etc. The order must be an even number 
!          and at least 4.
!       nfft1,nfft2,nfft3: the dimensions of the charge grid array
! OUTPUT
!       eer:  ewald reciprocal or k-space  energy
!       dx,dy,dz: forces incremented by k-space sum
!       virial:  virial due to k-space sum (valid for atomic scaling;
!                rigid molecule virial needs a correction term not
!                computed here
!-------------------------------------------------------------------------------------------

  integer numatoms,order,nfft1,nfft2,nfft3
  REAL*8 x(numatoms),y(numatoms),z(numatoms),charge(numatoms),charge0(numatoms),recip(3,3),volume, &
       ewald_coeff,rkcut

  REAL*8 eer,dx(numatoms),dy(numatoms),dz(numatoms),virial(3,3),dwrk

#ifndef FFTW
! SIZES of some arrays
  integer   sizfftab,sizffwrk,siztheta,siz_Q


! HEAP STORAGE:  These arrays need to be preserved throughout simulation
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3),fftable(sizfftab)
#endif
#ifdef FFTW 
! HEAP STORAGE:These arrays need to be preserved throughout simulation
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
! this can be tossed
  REAL*8, allocatable ::  Q(:,:,:),Q0(:,:,:)
  COMPLEX*16, allocatable :: FQ(:,:,:),FQ0(:,:,:)
  INTEGER ierr
  REAL*8 theta1(numatoms*order), theta2(numatoms*order),theta3(numatoms*order),dtheta1(numatoms*order), &
       dtheta2(numatoms*order),dtheta3(numatoms*order), fr1(numatoms), fr2(numatoms),fr3(numatoms),eer0
  type(C_PTR) planf,planr
#else
! STACK STORAGE: These arrays can be tossed after leaving this routine
  REAL*8 Q(siz_Q),Q0(siz_Q),ffwork(sizffwrk),theta1(siztheta),theta2(siztheta),theta3(siztheta), &
       dtheta1(siztheta),dtheta2(siztheta),dtheta3(siztheta),fr1(numatoms),fr2(numatoms),fr3(numatoms),eer0
  integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
!  get some integer array dimensions
#endif
#ifdef FFTW
  allocate(Q(nfft1,nfft2,nfft3),stat=ierr) 
  allocate(Q0(nfft1,nfft2,nfft3),stat=ierr) 
  allocate(FQ(1+nfft1/2,nfft2,nfft3),stat=ierr) 
  allocate(FQ0(1+nfft1/2,nfft2,nfft3),stat=ierr) 
  call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
  call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
! fill both l l0 gridded charge arrays 
  call alchemy_charge_grid(numatoms,charge,charge0,theta1,theta2,theta3,fr1,fr2,fr3,order, & 
       nfft1,nfft2,nfft3,Q,Q0)
!FT of gridded charge array at previous and current l point
  call dfftw_execute_dft_r2c(planf, Q, FQ)
  call dfftw_execute_dft_r2c(planf, Q0, FQ0)
! computes eer and eer0 
  call scalar_sum_alchemy(FQ,FQ0,ewald_coeff,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3, & 
       eer,eer0,virial,rkcut)
  dwrk = eer-eer0
! computes force normally 
  call dfftw_execute_dft_c2r(planr, FQ, Q)
  call grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,fr1,fr2,fr3, &
       order,nfft1,nfft2,nfft3,Q)
  deallocate(Q,Q0,FQ,FQ0,stat=ierr)
#else
  call get_fftdims(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
  call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
  call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
! fill both l l0 gridded charge arrays 
  call alchemy_charge_grid(numatoms,charge,charge0,theta1,theta2,theta3,fr1,fr2,fr3,order, & 
       nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,Q0)
!FT of gridded charge array at previous l point
  call fft_back(Q0,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
!FT of gridded charge array at current l point
  call fft_back(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
  call scalar_sum_alchemy(Q,Q0,ewald_coeff,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3, &
       nfftdim1,nfftdim2,nfftdim3,eer,eer0,virial,rkcut)
  dwrk = eer-eer0

!---  computes forces normally 
call fft_forward(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
call grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,fr1,fr2,fr3, &
     order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
#endif
return
end subroutine do_pmesh_alchemy

