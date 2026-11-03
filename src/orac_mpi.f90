module orac_mpi

! This version (NO_SYS_CALLS) compiles and runs on gfor 4.4 (lx01),
! Intel 10.1 (cl21), gfor 4.4 (fermi), bgxlf (fermi)
! added pre_dir
  
  use rem

#ifdef NO_SYS_CALLS
  use iso_c_binding
#endif
#ifdef _OMP_
#ifdef _BGQ_
      use omp_lib !omp layer  for BGQ only 
#endif
#endif


  implicit none
  include 'mpif.h'
  integer            :: iproc ! iproc is the MPI process rank (0-nproc-1)
  integer            :: nproc ! nproc is the total number of MPI processes
  integer            :: iproc_t !this PAR index of the process (1-nproc)
  integer            :: ntrajectories ! this is the number of trajectories (i.e. of PARxxx dirs)
  integer            :: nfile ! this is the number of systems to process

#ifdef NO_SYS_CALLS
  ! interface to standard C library functions 
  interface
     function mkdir(path, mode) bind(c)
       use iso_c_binding
       integer(c_int) mkdir
       character(kind=c_char) path
       integer(c_int16_t), value :: mode
     end function mkdir

     function getcwd(buffer,size) bind(c) result(r)
       use iso_c_binding
       character(kind=c_char) ,intent(out) :: buffer(*)
       integer(c_size_t),value,intent(in)  :: size
       type(c_ptr)  :: r
     end function getcwd

     function chdir(path) bind(c)
       use iso_c_binding
       integer(c_int) chdir
       character(kind=C_CHAR) path
     end function chdir

  end interface
#endif

contains 
  
  !-----------------------------------------------------------
  
  subroutine mpi_initialize
    !   initialize the MPI environment 
    implicit none
!    integer, intent(out) :: iproc
!    integer, intent(out) :: nproc
    ! local variables
    integer :: ierr
    
    CALL MPI_Init(ierr)
    
    CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    
    CALL MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)

    ! Before doing anything check whether nproc is less than maxrep
    if(nproc.gt.maxrep) THEN 
       if(iproc.eq.0) THEN
          WRITE(*,101) nproc, maxrep 
101       FORMAT( '************************************************' / &
               '*          FATAL ERROR!!!                      *' / &
               '*  maxrep is less than nproc                   *' / &
               '* maxstep =',i4,'  maxrep =',i4,'              *' / &
               '* Action: Increase maxrep in file rem.f90      *' / &
               '*         and recompile                        *' / &
               '************************************************' /)
          
          CALL MPI_Finalize(ierr)
          STOP
       end if
    end if
    nfile=0           ! initialize system count
  end subroutine mpi_initialize

!-----------------------------------------------------------


  subroutine mpi_mkdirs
    ! make a different dir for each trajectories identified by the color iproc_t
    ! and put the current wdir in it 
    implicit none
!    These variables  are globally defined in the module
!    integer,intent(in) :: iproc,nproc     ! this are ranks and size of MPI_COMM_WORLD
!    integer,intent(out) :: ntrajectories,iproc_t  ! this is the size and rank of trajectories (i.e. PARxxx dirs)

    ! local variables
    integer            :: ierr, scan, inull,nword,iret,clength
    character(len=128) :: directory,file,wdir,inputfile,pre_dir
    CHARACTER(len=130) :: line
    CHARACTER(len=1)   :: sep(2),comm(2)
    CHARACTER(len=80)  :: strngs(40),errmsg
    CHARACTER(len=9)   :: fmt
    LOGICAL            :: EXIST

#ifdef NO_SYS_CALLS
    integer(c_int)            :: i
    integer(c_int16_t) :: mode = o'0777'
    integer(c_size_t) :: size = 80
    type(c_ptr) buffer
#else
    integer            :: i
    integer            :: chdir, getcwd
    character(len=512) :: command
#endif
    sep(1)=" "
    sep(2)=","
    comm(1)="("
    comm(2)=")"

    !   copy the input file to a generic file INPUT (only rank0)
    !
    !   only rank0 process has access to _STDIN_, so it must first be
    !   copied to file INPUT, then each process must copy INPUT to
    !   PARXXXX/INPUT and read from there
    pre_dir="PAR"
    IF(iproc .EQ. 0) THEN
       OPEN (unit=99,file="INPUT")
       ntrajectories=nproc 
       rw_loop: do 
          READ(5,'(a130)',iostat=ierr) line
          if (line(1:3).EQ."#!&") THEN
             pre_dir=line(4:6)
          end if
          if(ierr < 0) exit rw_loop
          CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
          WRITE(99,'(a)') trim(line) !copy input in all cases 
       end do rw_loop
       CLOSE(99)
    END IF

! broadcast traj number to all processes       
    call MPI_BCAST(ntrajectories, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

! broadcast dirname to all processes       
    clength=len(pre_dir)
    call MPI_BCAST(pre_dir, clength, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr) 

    CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 
    
    ! initialize logical variables
    para_restart = .false.
    
    ! place the current working directory name in wdir
    wdir=''
    
#ifdef NO_SYS_CALLS
    buffer=getcwd(wdir,size)
    
#else
    ! getcwd is not a standard FORTRAN intrinsic 
    ! it is a language extension in gfortran, Intel
    ! i = 0 on success
    i=getcwd(wdir)

#endif


    ! if wdir ends with a NULL, remove it; 
    ! then, trim trailing blanks
    inull=SCAN(wdir, char(0), .true.) 
    if(inull.ne.0) wdir=wdir(1:inull-1)

    iproc_t=1+iproc  !assign replica/trajectory index (the iproc in the old code).

    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 
        ! directory='cwd/PARXXXX' ; Here XXXX is the number of trajectories 
    file=pre_dir
    WRITE (directory,'(a,a,a,i4.4)') trim(wdir),'/',trim(file),iproc_t

    ! check if the INPUT file is present in the PARXXXX directory
    exist=.false.
    inputfile = TRIM(directory) // '/INPUT'
    INQUIRE(FILE=inputfile,EXIST=exist) 

    ! try to change the work directory to PARXXXX
    ! i = 0 on success
    !   chdir is not a standard FORTRAN intrinsic 
    !   it is a language extension in gfortran, Intel
    i=chdir(trim(directory)//char(0))
    ! consistency check: 
    ! if there is a PARXXXX directory (i = 0), it must contain 
    ! an INPUT file from some old run
    if(i.eq.0.neqv.exist) THEN 
       write(6,*) ' Error: directory found, but no INPUT'
       write(6,*) ' process = ',iproc
       write(6,*) ' directory = <',trim(directory),'>, CHDIR return code=',i
       write(6,*) ' inputfile = <',inputfile,'> exist=',exist
       STOP
    END IF


    ! ***********************************************************************************
    ! ****************** if there is not a directory PARXXXX, create it and chdir to it *
    ! ***********************************************************************************
    if (i.ne.0) then 

       ! create directory PARXXXX

#ifdef NO_SYS_CALLS
       i=mkdir(trim(directory)//char(0),mode)
#else
       command = 'mkdir ' // directory // ' > /dev/null 2>&1'
       CALL system(trim(command))
#endif
       ! change the current work directory (anybody must do that) 
       ! i = 0 on success
       !   chdir is not a standard FORTRAN intrinsic 
       !   it is a language extension in gfortran, Intel.
       i=chdir(trim(directory)//char(0))
       if (i.ne.0) then
          ! error: stop
          WRITE(*,*) "************************************************"
          WRITE(*,*) "*          FATAL ERROR!!!                      *"
          WRITE(*,*) "*  Fatal error in changing directory           *"
          WRITE(*,*) "*  Check that the subdirectories PARXXXX       *"
          WRITE(*,*) "*  were created.                               *"
          WRITE(*,*) "************************************************"
          CALL MPI_Finalize(ierr)
          STOP 
       endif
100    continue
    endif

    ! ****************************************************************************************
    ! ************** the work directory is now PARXXXX                          **************
    ! ****************************************************************************************
    !   copy ../INPUT to INPUT
    OPEN (unit=98,file="../INPUT")
    OPEN (unit=99,file="INPUT")
    rw_loop_2: do 
       READ(98,'(a130)',iostat=ierr) line
       if(ierr < 0) exit rw_loop_2
       WRITE(99,'(a)') trim(line)
    end do rw_loop_2
    CLOSE(98)
    CLOSE(99)


    CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 

    ! refresh wdir
#ifdef NO_SYS_CALLS
    buffer=getcwd(wdir,size)
#else
    ! getcwd is not a standard FORTRAN intrinsic 
    ! it is a language extension in gfortran, Intel
    ! i = 0 on success
    i=getcwd(wdir)
#endif
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 
  end subroutine mpi_mkdirs

!-----------------------------------------------------------

  subroutine mpi_headings
#ifdef _OMP_
    use omp_integr, only:nthr
#endif
    implicit none
#ifdef _OMP_
#ifndef _BGQ_
    include 'omp_lib.h'  ! OMP layer
#endif
#endif
!   local variable 


#ifdef _OMP_
    if(iproc.eq.0) THEN
       write(*,123) ntrajectories,nthr
123    format &
            (10x,'========================================================' &
            /10x,'=                                                      =' &
            /10x,'=  O M P / M P I     O R A C     S T A R T E D         =' &
            /10x,'=                                                      =' &
            /10x,'=        TRAJECTORIES ARE ',i4,'                         =' &
            /10x,'=        FORCE PARALLELISM ON ',i4,' PROCESSORS          =' &
            /10x,'========================================================')
    end if
#else
    if(iproc.eq.0) THEN
       write(*,124) ntrajectories
124    format &
            (10x,'========================================================' &
            /10x,'=                                                      =' &
            /10x,'=          M P I     O R A C     S T A R T E D         =' &
            /10x,'=                                                      =' &
            /10x,'=        TRAJECTORIES ARE ',i4,'                         =' &
            /10x,'========================================================')
    end if
#endif
  end subroutine mpi_headings

!-----------------------------------------------------------

  subroutine mpi_close(iproc)
    implicit none
    integer, intent(in) :: iproc
    ! local variables 
    integer :: ierr
    
    if(iproc == 0) then 
       WRITE(6,100)
       WRITE(6,200)
    end if
    call MPI_FINALIZE(ierr)
    
100 FORMAT(// &
         '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/ &
         '     !                                                                      !'/ &
         '     !                        PROGRAM COMPLETED                             !'/ &
         '     !                                                                      !'/ &
         '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/) 
    
200 FORMAT(//10x, &
         '|--------------------------------------------------------------|',/10x, &
         "|    ``They're unfriendly.  Which is fortunate really;         |",/10x, &
         "|      they'd be difficult to love.''                          |",/10x, &
         "|                              . . . Avon    Blake's Seven     |",/10x, &
         "|                                                              |",/10x, &
         "|--------------------------------------------------------------|"//)
    
  end subroutine mpi_close

!-----------------------------------------------------------
  
end module orac_mpi
