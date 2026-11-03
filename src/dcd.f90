module dcd
  ! most of the code from lammps
  implicit none
  integer :: kdcd ! unit
  real(8) :: fdcd ! stride
  integer :: ndcd ! steps
  integer :: dcd_frame ! counter
  logical :: dcd_noh
  integer,allocatable :: list(:)
  
contains

  subroutine print_dcd(nato,box,x,y,z,beta)
    implicit none
    integer,intent(in) :: nato
    real(8),intent(in) :: box(3,3)
    real(8),intent(in) :: x(nato)
    real(8),intent(in) :: y(nato)
    real(8),intent(in) :: z(nato)
    character(7),intent(in) :: beta(nato)
    ! local variables
    integer :: j,k
    integer :: err
    real(8) :: a,b,c,aa,bb,gg
    integer :: nato1
    
    dcd_frame = dcd_frame + 1
    nato1 = nato
    
    if(dcd_frame == 1) then 
       if(dcd_noh) then 
          ! make a list of non-hydrogen atoms
          nato1 = count(beta(:)(1:1) /= 'h')
          allocate(list(nato1),stat=err)
          k = 0
          list = 0
          do j = 1,nato
             if(beta(j)(1:1) /= 'h') then 
                k = k + 1
                list(k) = j
             end if
          end do
       end if
       call print_dcd_header(nato1)
    end if

    a = box(1,1)
    c = sqrt(box(1,3)**2 + box(2,3)**2 + box(3,3)**2)
    bb = box(1,3) / c
    b = sqrt(box(1,2)**2 + box(2,2)**2)
    gg = box(1,2) / b
    aa = gg * (bb + box(2,3)*box(1,2)/ (c * box(2,2)))
    aa = acos(aa) / dacos(-1.d0) * 180.d0
    bb = acos(bb) / dacos(-1.d0) * 180.d0
    gg = acos(gg) / dacos(-1.d0) * 180.d0

    write(kdcd) 2*a, 2*b, 2*c, aa, bb, gg

    if(dcd_noh) then 
       write(kdcd) real(x((/ list /)))
       write(kdcd) real(y((/ list /)))
       write(kdcd) real(z((/ list /)))
    else
       write(kdcd) real(x)
       write(kdcd) real(y)
       write(kdcd) real(z)
    end if
    
  end subroutine print_dcd
  
  subroutine print_dcd_header(nato)
    implicit none
    integer,intent(in) :: nato
    character*4 :: hdr
    integer :: icntrl(20), nstr
    integer,parameter :: max_dumps=1000
    
    hdr = 'CORD'
    icntrl = 0
    nstr = 0! number of strings in header
    icntrl(1) = max_dumps! number of frames in traj file
    icntrl(2) = 0! number of steps in previous run
    icntrl(3) = 1! frequency of saving
    icntrl(4) = max_dumps! total number of steps
    icntrl(8) = nato*3 - 6! number of degrees of freedom
    icntrl(10) = 981668463! coded time step
    icntrl(11) = 1! coded crystallographic group (or zero)
    icntrl(20) = 28! CHARMM version number
    
    write(kdcd) hdr, icntrl
    write(kdcd) nstr
    write(kdcd) nato
    
  end subroutine print_dcd_header

end module dcd
