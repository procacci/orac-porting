#ifndef FFTW
subroutine pmesh_kspace_get_sizes(nfft1,nfft2,nfft3,numatoms,order,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
  implicit none
  integer nfft1,nfft2,nfft3,numatoms,order, sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack

  ! INPUT
  !nfft1,nfft2,nfft3,numatoms,order
  !nfft1,nfft2,nfft3 are the dimensions of the charge grid array
  !numatoms is number of atoms
  !order is the order of B-spline interpolation

  ! OUTPUT
  !sizfftab,sizffwrk,siztheta,siz_Q
  !sizfftab is permanent 3d fft table storage
  !sizffwrk is temporary 3d fft work storage
  !siztheta is size of arrays theta1-3 dtheta1-3
  !sizheap is total size of permanent storage
  !sizstack is total size of temporary storage


  ! This routine computes the above output parameters needed for 
  ! heap or stack allocation.

  integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork

  call get_fftdims(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sizfftab,sizffwrk)
  siztheta = numatoms*order
  siz_Q = 2*nfftdim1*nfftdim2*nfftdim3
  sizheap = nfft1+nfft2+nfft3+sizfftab
  sizstack = siz_Q+6*siztheta+sizffwrk+3*numatoms
  !write(6,*)'total HEAP storage needed = ',sizheap
  !write(6,*)'total STACK storage needed = ',sizstack
  return
end subroutine pmesh_kspace_get_sizes
!----------------------------------------------------
subroutine pmesh_kspace_setup( bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork, nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
  implicit none

  !see DO_PMESH_KSPACE for explanation of arguments

  integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
  REAL*8 fftable(sizfftab),ffwork(sizffwrk)

  REAL*8 dummy
  integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

  call get_fftdims(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
  call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,order)
  call fft_setup(dummy,fftable,ffwork, nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, nfftable,nffwork)
  return
end subroutine pmesh_kspace_setup
!----------------------------------------------------
#endif
#ifdef FFTW 
subroutine do_pmesh_kspace(numatoms,x,y,z,charge,recip,volume,ewald_coeff,order,nfft1,nfft2,nfft3,eer,dx,dy,dz,virial,bsp_mod1,bsp_mod2,bsp_mod3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fr1,fr2,fr3,rkcut,planf,planr)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'  
#else
subroutine do_pmesh_kspace(numatoms,x,y,z,charge,recip,volume,ewald_coeff,order,nfft1,nfft2,nfft3,eer,dx,dy,dz,virial,sizfftab,sizffwrk,siztheta,siz_Q,bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fr1,fr2,fr3,rkcut)
  implicit none
#endif
  ! INPUT 
  ! numatoms:number of atoms
  ! x,y,z: atomic coords
  ! chargeatomic charges
  ! recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
  ! volume: the volume of the unit cell
  ! ewald_coeff: ewald convergence parameter
  ! order: the order of Bspline interpolation. E.g. cubic is order 4
  !fifth degree is order 6 etc. The order must be an even number 
  !and at least 4.
  ! nfft1,nfft2,nfft3: the dimensions of the charge grid array
  integer numatoms,order,nfft1,nfft2,nfft3
  REAL*8 x(numatoms),y(numatoms),z(numatoms),charge(numatoms),recip(3,3),volume,ewald_coeff,rkcut

  ! OUTPUT
  ! eer:ewald reciprocal or k-spaceenergy
  ! dx,dy,dz: forces incremented by k-space sum
  ! virial:virial due to k-space sum (valid for atomic scaling;
  !rigid molecule virial needs a correction term not
  !computed here
  REAL*8 eer,dx(numatoms),dy(numatoms),dz(numatoms), virial(3,3)

#ifndef FFTW
  ! SIZES of some arrays
  integer sizfftab,sizffwrk,siztheta,siz_Q

  ! HEAP STORAGE:These arrays need to be preserved throughout simulation
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3),fftable(sizfftab)
#endif
  ! STACK STORAGE: These arrays can be tossed after leaving this routine
#ifdef FFTW 
  ! this can be tossed
  real*8, allocatable ::  Q(:,:,:)
  complex*16, allocatable :: FQ(:,:,:)
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
  real*8 theta1(numatoms*order), theta2(numatoms*order),theta3(numatoms*order),dtheta1(numatoms*order), dtheta2(numatoms*order),dtheta3(numatoms*order),fr1(numatoms), fr2(numatoms),fr3(numatoms)
  integer ierr
  type(C_PTR) planf,planr
#else
  REAL*8 Q(siz_Q),ffwork(sizffwrk),theta1(siztheta), theta2(siztheta),theta3(siztheta),dtheta1(siztheta), dtheta2(siztheta),dtheta3(siztheta),fr1(numatoms), fr2(numatoms),fr3(numatoms)
  integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
  !get some integer array dimensions
#endif

#ifdef FFTW
  allocate(Q(nfft1,nfft2,nfft3),stat=ierr) 
  allocate(FQ(1+nfft1/2,nfft2,nfft3),stat=ierr) 
  call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
  call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
  call fill_charge_grid(numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
  call fftw_execute_dft_r2c(planf, Q, FQ)
  call scalar_sum(FQ,ewald_coeff,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,eer,virial,rkcut)
  call fftw_execute_dft_c2r(planr, FQ, Q)
  call grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
deallocate(Q,FQ,stat=ierr)
#else
  call get_fftdims(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
  call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
  call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
  call fill_charge_grid(numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
  call fft_back(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
  call scalar_sum(Q,ewald_coeff,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,virial,rkcut)
  call fft_forward(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
  call grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
#endif
  return
end subroutine do_pmesh_kspace
!----------------------------------------------------------------------
subroutine get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
  implicit none

  ! INPUT:
  !numatoms: number of atoms
  !x,y,z: arrays of cartesian coords
  !recip: the 3x3 array of reciprocal vectors stored as columns
  ! OUTPUT:
  ! fr1,fr2,fr3 the scaled and shifted fractional coords

  integer numatoms,nfft1,nfft2,nfft3
  REAL*8 x(numatoms),y(numatoms),z(numatoms),recip(3,3)
  REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)

  integer n
  REAL*8 w
#ifdef _OMP_
!$OMP  PARALLEL DO PRIVATE(w) SCHEDULE(STATIC) 
#endif
  do  n = 1,numatoms
     w = x(n)*recip(1,1)+y(n)*recip(2,1)+z(n)*recip(3,1)
     fr1(n) = nfft1*(w - anint(w) + 0.5d0)
     w = x(n)*recip(1,2)+y(n)*recip(2,2)+z(n)*recip(3,2)
     fr2(n) = nfft2*(w - anint(w) + 0.5d0)
     w = x(n)*recip(1,3)+y(n)*recip(2,3)+z(n)*recip(3,3)
     fr3(n) = nfft3*(w - anint(w) + 0.5d0)
  end do
#ifdef _OMP_
!$OMP END PARALLEL DO
#endif
  return
end subroutine get_scaled_fractionals
!---------------------------------------------------------------
subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,order)
  implicit none
  integer nfft1,nfft2,nfft3,order
  REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)

  integer MAXORDER
  parameter (MAXORDER=25)
  integer MAXNFFT
  parameter (MAXNFFT=1000)
  REAL*8 array(MAXORDER),darray(MAXORDER),w
  REAL*8 bsp_arr(MAXNFFT)
  integer i,maxn

  ! this routine loads the moduli of the inverse DFT of the B splines
  ! bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
  ! Order is the order of the B spline approx.

  if ( order .gt. MAXORDER )then
     write(6,*)'order too large! check on MAXORDER'
     stop
  endif
  maxn = max(nfft1,nfft2,nfft3)
  if ( maxn .gt. MAXNFFT )then 
     write(6,*)'nfft1-3 too large! check on MAXNFFT'
     stop
  endif
  w = 0.d0
  call fill_bspline(w,order,array,darray)
  do  i = 1,maxn
     bsp_arr(i) = 0.d0
  end do
  do  i = 2,order+1
     bsp_arr(i) = array(i-1)
  end do
  call DFTMOD(bsp_mod1,bsp_arr,nfft1)
  call DFTMOD(bsp_mod2,bsp_arr,nfft2)
  call DFTMOD(bsp_mod3,bsp_arr,nfft3)
  return
end subroutine load_bsp_moduli
!------------------------------------------------------------------------
subroutine DFTMOD(bsp_mod,bsp_arr,nfft)
  implicit none
  integer nfft
  REAL*8 bsp_mod(nfft),bsp_arr(nfft)
  ! Computes the modulus of the discrete fourier transform of bsp_arr,
  !storing it into bsp_mod

  integer j,k
  REAL*8 sum1,sum2,twopi,arg,tiny
  twopi = 2.d0*3.14159265358979323846
  tiny = 1.d-7
  do  k = 1,nfft
     sum1 = 0.d0
     sum2 = 0.d0
     do  j = 1,nfft
        arg = twopi*(k-1)*(j-1)/nfft
        sum1 = sum1 + bsp_arr(j)*dcos(arg)
        sum2 = sum2 + bsp_arr(j)*dsin(arg)
     end do
     bsp_mod(k) = sum1**2 + sum2**2
  end do
  do  k = 1,nfft
     if ( bsp_mod(k) .lt. tiny )bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
  end do
  return
end subroutine DFTMOD
!------------------------------------------------------------------------

#ifdef FFTW
subroutine fill_charge_grid(numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
#else
subroutine fill_charge_grid(numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
#endif
  !---------------------------------------------------------------------
  ! INPUT:
  !numatoms:number of atoms
  !charge: the array of atomic charges
  !theta1,theta2,theta3: the spline coeff arrays
  !fr1,fr2,fr3 the scaled and shifted fractional coords
  !nfft1,nfft2,nfft3: the charge grid dimensions
  !nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
  !order: the order of spline interpolation
  ! OUTPUT:
  !Q the charge grid
  !---------------------------------------------------------------------
  implicit none
  integer numatoms,order,nfft1,nfft2,nfft3
  integer nfftdim1,nfftdim2,nfftdim3
  REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)
  REAL*8 theta1(order,numatoms),theta2(order,numatoms),theta3(order,numatoms),charge(numatoms)
#ifdef FFTW
  REAL*8 Q(nfft1,nfft2,nfft3)
#else
  REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3)
#endif
  integer n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
  REAL*8 prod

  Q=0  ! clear Q whatever the dimensions are 

#ifdef _OMP_
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) & 
!$OMP& SHARED(numatoms,fr1,fr2,fr3,charge,order) &
!$OMP& SHARED(nfft1,nfft2,nfft3,theta1,theta2,theta3) &
!$OMP& REDUCTION(+:Q)  
#endif
  do n = 1,numatoms
     k0 = int(fr3(n)) - order
     do  ith3 = 1,order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = int(fr2(n)) - order
        do  ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod = theta2(ith2,n)*theta3(ith3,n)*charge(n)
           i0 = int(fr1(n)) - order
           do  ith1 = 1,order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
#ifdef FFTW               
              Q(i,j,k) = Q(i,j,k) + theta1(ith1,n) * prod
#else
              Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod
#endif
           end do
        end do
     end do
  end do
#ifdef _OMP_
!$OMP END PARALLEL DO 
#endif 
  return
end subroutine fill_charge_grid
!-----------------------------------------------------------
subroutine clearQ(Q,ntot)
  integer ntot
  REAL*8 Q(ntot)
  integer i
  do i = 1,ntot
     Q(i) = 0.d0
  end do
  return
end subroutine clearQ
!-----------------------------------------------------------
#ifdef FFTW
subroutine grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fx,fy,fz,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
#else
subroutine grad_sum(numatoms,charge,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fx,fy,fz,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
#endif
  implicit none
  integer numatoms,order,nfft1,nfft2,nfft3
  integer nfftdim1,nfftdim2,nfftdim3
  REAL*8 recip(3,3)
  REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)
  REAL*8 fx(numatoms),fy(numatoms),fz(numatoms)
  REAL*8 theta1(order,numatoms),theta2(order,numatoms),theta3(order,numatoms),charge(numatoms)
  REAL*8 dtheta1(order,numatoms),dtheta2(order,numatoms),dtheta3(order,numatoms)
#ifdef FFTW
  REAL*8 Q(nfft1,nfft2,nfft3)
#else
  REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3)
#endif
  integer n,ith1,ith2,ith3,i0,j0,k0,i,j,k
  REAL*8 f1,f2,f3,term

#ifdef _OMP_
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC)  & 
!$OMP& SHARED(numatoms,fr1,fr2,fr3,charge,Q,fx,fy,fz,recip,order) & 
!$OMP& SHARED(nfft1,nfft2,nfft3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
#endif
  do n = 1,numatoms
     f1 = 0.d0
     f2 = 0.d0
     f3 = 0.d0
     k0 = int(fr3(n)) - order
     do  ith3 = 1,order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = int(fr2(n)) - order
        do  ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           i0 = int(fr1(n)) - order
           do  ith1 = 1,order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
#ifdef FFTW
              term = charge(n)*Q(i,j,k)
#else
              term = charge(n)*Q(1,i,j,k)
#endif
              ! force is negative of grad
              f1 = f1 - nfft1 * term * dtheta1(ith1,n) * theta2(ith2,n) * theta3(ith3,n)
              f2 = f2 - nfft2 * term * theta1(ith1,n) * dtheta2(ith2,n) * theta3(ith3,n)
              f3 = f3 - nfft3 * term * theta1(ith1,n) * theta2(ith2,n) * dtheta3(ith3,n)
           end do
        end do
     end do
     fx(n) = fx(n) + recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
     fy(n) = fy(n) + recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
     fz(n) = fz(n) + recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
  end do
#ifdef _OMP_
!$OMP END PARALLEL DO 
#endif 
  return
end subroutine grad_sum

!-------------------------------------------------------------
subroutine check_virial(self_ene,adj_ene,dir_ene,rec_ene,adj_vir,rec_vir,dir_vir)
  implicit none
  REAL*8 self_ene,adj_ene,dir_ene,rec_ene
  REAL*8 adj_vir(6),rec_vir(6),dir_vir(6)

  REAL*8 etot,svir,relerr
  etot = self_ene+adj_ene+dir_ene+rec_ene
  svir = adj_vir(1)+rec_vir(1)+dir_vir(1)+adj_vir(4)+rec_vir(4)+dir_vir(4)+adj_vir(6)+rec_vir(6)+dir_vir(6)
  relerr = 2.d0*abs(etot+svir)/(abs(etot)+abs(svir))
  write(6,*)'tot ene = ',etot
  write(6,*)'trace vir = ',svir
  write(6,*)'rel error = ',relerr
  return
end subroutine check_virial
!-------------------------------------------------------------
subroutine check_force(numatoms, fx1,fy1,fz1,fx2,fy2,fz2,fdx,fdy,fdz)
  integer numatoms
  REAL*8 fx1(*),fy1(*),fz1(*),fx2(*),fy2(*),fz2(*), fdx(*),fdy(*),fdz(*)

  REAL*8 rms_num,rms_den,rms
  integer i
  rms_num = 0.d0
  rms_den = 0.d0
  do  i = 1,numatoms
     rms_num = rms_num + (fx1(i)-fx2(i))**2 + (fy1(i)-fy2(i))**2 + (fz1(i)-fz2(i))**2
     rms_den = rms_den + fdx(i)**2 + fdy(i)**2 + fdz(i)**2
  end do
  rms = dsqrt(rms_num/rms_den)
  write(6,*)'rms force err = ',rms
  return
end subroutine check_force
!-------------------------------------------------------------

subroutine pubz3di(n1,n2,n3,table,ntable)
  implicit none
  integer n1,n2,n3,ntable
  REAL*8 table(ntable,3)
  ! ntable should be 4*max(n1,n2,n3) +15


  call cffti(n1,table(1,1))
  call cffti(n2,table(1,2))
  call cffti(n3,table(1,3))

  return
end subroutine pubz3di
!*****************************************************************************
subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable, work,nwork)
  implicit none

  integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
  COMPLEX*16 w(ld1,ld2,n3)
  COMPLEX*16 work( nwork)
  REAL*8 table(ntable,3)

  integer i,j,k
  ! ntable should be 4*max(n1,n2,n3) +15
  ! nwork should be max(n1,n2,n3)

  ! transform along Xfirst ...

  do  k = 1, n3
     do  j = 1, n2
        do  i = 1,n1
           work(i) = w(i,j,k)
        end do
        if ( isign .eq. -1) call cfftf(n1,work,table(1,1))
        if ( isign .eq. 1) call cfftb(n1,work,table(1,1))
        do  i = 1,n1
           w(i,j,k) = work(i)
        end do
     end do
  end do

  ! transform along Y then ...

  do  k = 1,n3
     do  i = 1,n1
        do  j = 1,n2
           work(j) = w(i,j,k)
        end do
        if ( isign .eq. -1) call cfftf(n2,work,table(1,2))
        if ( isign .eq. 1) call cfftb(n2,work,table(1,2))
        do  j = 1,n2
           w(i,j,k) = work(j)
        end do
     end do
  end do

  ! transform along Z finally ...

  do  i = 1, n1
     do  j = 1, n2
        do  k = 1,n3
           work(k) = w(i,j,k)
        end do
        if ( isign .eq. -1) call cfftf(n3,work,table(1,3))
        if ( isign .eq. 1) call cfftb(n3,work,table(1,3))
        do  k = 1,n3
           w(i,j,k) = work(k)
        end do
     end do
  end do

  return
end subroutine pubz3d
!----------------------------------------------------


! SPLINE

!---------------------------------------------------------------------
subroutine get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
  !---------------------------------------------------------------------
  ! INPUT:
  !numatoms: number of atoms
  !fr1,fr2,fr3 the scaled and shifted fractional coords
  !order: the order of spline interpolation
  ! OUTPUT
  !theta1,theta2,theta3: the spline coeff arrays
  !dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
  !---------------------------------------------------------------------
  implicit none
  integer numatoms,order
  REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)
  REAL*8 theta1(order,numatoms),theta2(order,numatoms),theta3(order,numatoms),dtheta1(order,numatoms),dtheta2(order,numatoms),dtheta3(order,numatoms)

  REAL*8 w
  integer n

!!$OMP PARALLEL DO private(w) SCHEDULE(STATIC)
  do  n = 1,numatoms
     w = fr1(n)-int(fr1(n))
     call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
     w = fr2(n)-int(fr2(n))
     call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
     w = fr3(n)-int(fr3(n))
     call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
  end do
!!$OMP END PARALLEL DO
  return
end subroutine get_bspline_coeffs
!---------------------------------------------------
subroutine fill_bspline(w,order,array,darray)
  !---------- use standard B-spline recursions: see doc file
  implicit none
  integer order
  REAL*8 w,array(order),darray(order)

  integer k
  ! do linear case
  call init(array,w,order)
  ! compute standard b-spline recursion
  do k = 3,order-1
     call one_pass(array,w,k)
  end do
  ! perform standard b-spline differentiation
  call diff(array,darray,order)
  ! one more recursion
  call one_pass(array,w,order)
  return
end subroutine fill_bspline
!---------------------------------------------------
subroutine init(c,x,order)
  implicit none
  integer order
  REAL*8 c(order),x
  c(order) = 0.d0
  c(2) = x
  c(1) = 1.d0 - x
  return
end subroutine init
!-------------------------------------
subroutine one_pass(c,x,k)
  implicit none
  REAL*8 c(*),x
  integer k

  REAL*8 div
  integer j

  div = 1.d0 / (k-1)
  c(k) = div*x*c(k-1)
  do  j = 1,k-2
     c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
  end do
  c(1) = div*(1-x)*c(1)
  return
end subroutine one_pass
!-------------------------------------
subroutine diff(c,d,order)
  implicit none
  REAL*8 c(*),d(*)
  integer order

  integer j
  d(1) = -c(1)
  do  j = 2,order
     d(j) = c(j-1) - c(j)
  end do
  return
end subroutine diff
!-------------------------------------


! FFT CALLS

subroutine get_fftdims(nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sizfftab,sizffwrk)
  implicit none
  integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sizfftab,sizffwrk
  integer n,nfftmax

#ifdef FTTW
  nfftdim1=nfft1 
  nfftdim2=nfft2 
  nfftdim3=1+nfft3/2   ! this is done for in-place t2c or c2r FFT
#else
  nfftmax = max(nfft1,nfft2,nfft3)
  nfftdim1 = nfft1
  n = nfft1/2
  if ( nfft1 .eq. 2*n )nfftdim1 = nfft1+1
  nfftdim2 = nfft2
  n = nfft2/2
  if ( nfft2 .eq. 2*n )nfftdim2 = nfft2+1
  nfftdim3 = nfft3
  n = nfft3/2
  if ( nfft3 .eq. 2*n )nfftdim3 = nfft3+1
  nfftable = 4*nfftmax + 15
  nffwork = nfftmax
  sizfftab = 3*nfftable
  sizffwrk= 2*nfftmax
#endif

  return
end subroutine get_fftdims

subroutine fft_setup(array,fftable,ffwork, nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, nfftable,nffwork)
  implicit none

  REAL*8 array(*),fftable(*),ffwork(*)
  integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  integer nfftable,nffwork,isys(4)

  integer isign,inc1,inc2,inc3
  REAL*8 scale

  call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)

  return
end subroutine fft_setup
!-----------------------------------------------------------
subroutine fft_forward(array,fftable,ffwork, nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, nfftable,nffwork)
  implicit none

  REAL*8 array(*),fftable(*),ffwork(*)
  integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

  integer isign,inc1,inc2,inc3
  REAL*8 scale
  integer nfftable,nffwork,isys(4)

  isign = 1

  call pubz3d(isign,nfft1,nfft2,nfft3,array,nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

  return
end subroutine fft_forward
!-----------------------------------------------------------
subroutine fft_back(array,fftable,ffwork, nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, nfftable,nffwork)
  implicit none

  REAL*8 array(*),fftable(*),ffwork(*)
  integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  integer nfftable,nffwork,isys(4)

  integer isign,inc1,inc2,inc3
  REAL*8 scale

  isign = -1

  call pubz3d(isign,nfft1,nfft2,nfft3,array,nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

  return
end subroutine fft_back

