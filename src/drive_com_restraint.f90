!---------------------------------------------------------------------------------
subroutine com_restraint (k,r0,x,y,z,m1,m2,m,r,fx,fy,fz,u,ilig1,ilig2,itar1,itar2,lforce)

  implicit none

! Variables included in the call
  real(8) :: k, r0, u, r 
  real(8) :: x(*), y(*), z(*), m(*), m1, m2, fx(*), fy(*), fz(*)
  integer(4) :: ilig1,ilig2,itar1,itar2
  LOGICAL :: lforce

! Local variables
  real(8) :: qforce,x1,y1,z1,x2,y2,z2,r1,x21,y21,z21
  real(8) :: fx1,fy1,fz1,fx2,fy2,fz2
  integer :: i

! Coordinates of the centers of mass of the fragments
  x1 = 0.
  y1 = 0.
  z1 = 0.
  x2 = 0.
  y2 = 0.
  z2 = 0.
  do i = ilig1, ilig2
     x1 = x1 + m(i) * x(i)
     y1 = y1 + m(i) * y(i)
     z1 = z1 + m(i) * z(i)
  end do
!$OMP  PARALLEL DO 
  do i = itar1, itar2
     x2 = x2 + m(i) * x(i)
     y2 = y2 + m(i) * y(i)
     z2 = z2 + m(i) * z(i)
  end do
!$OMP END PARALLEL DO 
  x1 = x1*m1
  y1 = y1*m1
  z1 = z1*m1
  x2 = x2*m2
  y2 = y2*m2
  z2 = z2*m2

! Distance between the centers of mass of the fragments
  r = sqrt( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2 )
  IF(.not.lforce) RETURN
  r1=1.d0/r
  x21=x2-x1
  y21=y2-y1
  z21=z2-z1

! Calculation of the force on the centers of mass (divided by the frags masses)
  qforce=2.0D0*k*(r-r0)*r1

  fx1 = qforce*x21*m1
  fy1 = qforce*y21*m1
  fz1 = qforce*z21*m1

  fx2 = -qforce*x21*m2
  fy2 = -qforce*y21*m2
  fz2 = -qforce*z21*m2


! Distribution of the forces on the atoms of the fragments
  do i = ilig1,ilig2
     fx(i) = fx(i) + m(i) * fx1
     fy(i) = fy(i) + m(i) * fy1
     fz(i) = fz(i) + m(i) * fz1
  end do
!$OMP  PARALLEL DO 
  do i = itar1,itar2
     fx(i) = fx(i) + m(i) * fx2
     fy(i) = fy(i) + m(i) * fy2
     fz(i) = fz(i) + m(i) * fz2
  enddo
!$OMP END PARALLEL DO 

! Potential energy
  u = k * ( r - r0 )**2
  return
end subroutine com_restraint

subroutine comp_work_com(rcom,r0com,kcom,str_com,t,t0,t1,tstep,workcom,diff1,nplotstn0,ninn0,kplotst0,iaux) 

  implicit NONE
  real*8 rcom,r0com,kcom,workcom,str_com,tstep,t,t0,t1 ! arguments
  integer nplotstn0,ninn0,kplotst0,iaux                ! arguments
  real*8 diff0,diff1,xkv,lambda
  logical steering
  if(t.lt.t0.or.t.gt.t1) RETURN
  diff0=diff1
  xkv = kcom*str_com*4.184 
  diff1= (rcom-r0com)
  !In the funcion "work" there is a 2 factor from the 1st der. that
  !is canceling with a 1/2 factor from trapez. numeric. integration 
  workcom = workcom - xkv*(diff1+diff0)*tstep
  IF(mod(ninn0,nplotstn0).EQ.0) THEN 
     if(iaux.lt.0) THEN 
        write(kplotst0,1015) t,r0com,rcom,workcom
     else
        write(kplotst0,1016) t,r0com,rcom,workcom,iaux
     endif
1015 FORMAT(" COM ", g15.3,3g15.5) 
1016 FORMAT(" COM ", g15.3,3g15.5, "  Replica -> ",i4)   
  END IF
  
end subroutine comp_work_com
