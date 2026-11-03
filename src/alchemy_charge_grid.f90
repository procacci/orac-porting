#ifdef FFTW
subroutine alchemy_charge_grid(numatoms,charge,charge0,theta1,theta2,theta3,fr1,fr2,fr3, &
     order,nfft1,nfft2,nfft3,Q,Q0)
#else
subroutine alchemy_charge_grid(numatoms,charge,charge0,theta1,theta2,theta3,fr1,fr2,fr3, &
     order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,Q0)
#endif
    !---------------------------------------------------------------------
! INPUT:
!      numatoms:  number of atoms
!      charge: the array of atomic charges
!      theta1,theta2,theta3: the spline coeff arrays
!      fr1,fr2,fr3 the scaled and shifted fractional coords
!      nfft1,nfft2,nfft3: the charge grid dimensions
!      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
!      order: the order of spline interpolation
! OUTPUT:
!      Q the charge grid
!                    !!!!! WARNING  !!!!!!
! with intel compile for huge arrays stack can give problems with OpenMP
! Use export KMP_STACKSIZE=1g before compiling
!                    !!!!! WARNING  !!!!!!
!
!---------------------------------------------------------------------
  implicit none
  integer numatoms,order,nfft1,nfft2,nfft3
  REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)
  REAL*8 theta1(order,numatoms),theta2(order,numatoms),theta3(order,numatoms), & 
       charge(numatoms),charge0(numatoms)
#ifdef FFTW
  real*8 Q(nfft1,nfft2,nfft3),Q0(nfft1,nfft2,nfft3),prod0
#else
  REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3),Q0(2,nfftdim1,nfftdim2,nfftdim3),prod0
#ifdef _OMP_
  REAL*8, allocatable :: Q0T(:,:,:),QT(:,:,:)
#endif
  integer nfftdim1,nfftdim2,nfftdim3,ierr
#endif
  integer n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
  REAL*8 prod
  Q=0  ! clear Q 
  Q0=0  ! clear Q0 
#ifndef FFTW
#ifdef _OMP_
if(.not.allocated(QT)) allocate(QT(nfft1,nfft2,nfft3),stat=ierr)
if(.not.allocated(Q0T)) allocate(Q0T(nfft1,nfft2,nfft3),stat=ierr)
  QT=0  ! clear Q 
  Q0T=0  ! clear Q0 
#endif
#endif
  
#ifdef _OMP_
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) & 
!$OMP& SHARED(numatoms,fr1,fr2,fr3,charge,charge0,order) &
!$OMP& SHARED(nfft1,nfft2,nfft3,theta1,theta2,theta3) &
#ifndef FFTW
!$OMP& REDUCTION(+:QT,Q0T)  
#else
!$OMP& REDUCTION(+:Q,Q0)  
#endif
#endif
  do n = 1,numatoms
     k0 = int(fr3(n)) - order
#ifndef _BGQ_
#ifndef _PGI_
     if(isnan(fr3(n))) THEN 
        write(6,*) "fr3 in alchemy charge_grid is a NaN",n
        STOP
     end if
#endif
#endif
     do ith3 = 1,order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        j0 = int(fr2(n)) - order
        do ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           prod = theta2(ith2,n)*theta3(ith3,n)*charge(n)
           prod0 = theta2(ith2,n)*theta3(ith3,n)*charge0(n)
           i0 = int(fr1(n)) - order
           do ith1 = 1,order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
#ifdef FFTW
              Q(i,j,k) = Q(i,j,k) + theta1(ith1,n) * prod
              Q0(i,j,k) = Q0(i,j,k) + theta1(ith1,n) * prod0
#else
#ifdef _OMP_
              QT(i,j,k) = QT(i,j,k) + theta1(ith1,n) * prod
              Q0T(i,j,k) = Q0T(i,j,k) + theta1(ith1,n) * prod0
#else
              Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod
              Q0(1,i,j,k) = Q0(1,i,j,k) + theta1(ith1,n) * prod0
#endif
#endif
           END DO
        END DO
     END DO
  END DO
#ifdef _OMP_
!$OMP END PARALLEL DO 
#ifndef FFTW
!$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(Q,QT,Q0,Q0T,nfft1,nfft2,nfft3)
  do k=1,nfft3
     do j=1,nfft2
        do i=1,nfft1
           Q(1,i,j,k)=QT(i,j,k) 
           Q(2,i,j,k)=0.d0
           Q0(1,i,j,k)=Q0T(i,j,k) 
           Q0(2,i,j,k)=0.d0
        end do
     end do
  end do
!$OMP END PARALLEL DO
#endif
#endif 
  return
end subroutine alchemy_charge_grid
