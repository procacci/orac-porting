module rem 
  use dcd, only: ndcd
  implicit none
  INTEGER, parameter :: maxrep=1024,maxseg=30 ! max values of replica number and segment number 
  INTEGER, parameter :: potterms=3           ! number of terms in total potential subdivision	 
  INTEGER :: rem_group(2,maxseg)             ! matrix of the first and the last atoms of each group
  INTEGER :: rem_acc(maxrep)                 ! number of accepted exchanges in replica index order
  INTEGER :: rep_acc(maxrep)                 ! number of accepted exchanges in replica index order (for the current process)
  INTEGER :: rem_table(0:maxrep-1)           ! array of replica indeces in processor index order 
  INTEGER :: para_index        ! processor index
  INTEGER :: rem_groups        ! total number of REM segments
  INTEGER :: rem_segkind       ! 1: only the segments-segments interactions are divided; 2: segments-system interactions; 0: both.
  INTEGER :: cdist             ! input parameter. 0: restart; 1 default temperature progression; 2 user defined temperature progression
  INTEGER :: step_rem          ! total REM exchanges counter 
  INTEGER :: krem              ! output unit for rem energies
  integer :: nrem
  integer :: krem_dcd
  integer :: nbatteries        ! number of identical batteries (same scaling x nbatteries time) 
  real*8  :: frem              ! frequency for rem data printing
  REAL*8  :: rem_mat(maxrep,potterms) ! matrix of scaling factors ordered for replica index
  REAL*8  :: rem_factor(potterms)     ! array of multiplicative scaling factor
  REAL*8  :: rem_pot(potterms)        ! array of (unscaled) potential energies (n0,n1,nb)
  REAL*8  :: rem_factor_max(potterms) ! array of maximum values of the factors
  REAL*8  :: rem_print         ! print interval
  REAL*8  :: rem_adjust_ratio  ! maximum simulation time for exchange ratio adjustment
  REAL*8  :: rem_printd        ! print interval REM_DIAGNOSTIC
  REAL*8  :: rem_ts            ! exchange interval
  REAL*8  :: rem_n0_energy     ! (unscaled) i-torsional energy of the selected segments
  REAL*8  :: rem_n1_energy     ! (unscaled) torsional energy of the selected segments
  REAL*8  :: rem_hnb_energy    ! (unscaled) h-shell energy of the selected segments
  REAL*8  :: rem_lnb_energy    ! (unscaled) l-shell energy of the selected segments
  REAL*8  :: rem_mnb_energy    ! (unscaled) m-shell energy of the selected segments
  LOGICAL :: para_restart      ! if true, restart from a previous parallel run
  LOGICAL :: rem_run           ! if true, a REMD run is performed
  LOGICAL :: rem_segment       ! if true, only the torsional potential of the selected segments is scaled
  LOGICAL :: prnt_remdt  ! if true, print rem energies with the same frequency of pdb
  LOGICAL :: lskip_force  ! if true, force routine is skipped. This is useful when sampling in the
                         ! gas phase with INTRAMOL=T.  
 
  character(80) :: dcd_rem_file

contains 

!-----------------------------------------------------------

  subroutine print_remdata(dt,ninner,mrespa,lrespa,para0)
    implicit none
    integer mod
    real(8),intent(in) :: dt
    integer, intent(in) :: ninner,mrespa,lrespa,para0
    real(8) :: fstep

    fstep=dt*DFLOAT(ninner)/dfloat(mrespa*lrespa)
    ! print rem data 
    if(prnt_remdt) then 
       if(nrem /= 0) then 
          if(mod(ninner,nrem) == 0) then 
             write(krem,'(f15.1,i10,3f15.5,i10)') fstep,para_index,rem_pot,para0
          end if
       end if
    end if
    
    ! print rem data whit the same frequency of dcd trajectory
    if(ndcd /= 0) then 
       if(mod(ninner,ndcd) == 0) then 
          write(krem_dcd,'(f15.1,i10,3f15.5,i10)') fstep,para_index,rem_pot,para0
       end if
    end if

  end subroutine print_remdata

#ifdef _MPI_

!-----------------------------------------------------------
  
  ! todo: move remltor,remint14,rematom in the module
  subroutine init_rem(iproc,nproc,kprint,ntap,ltors,litor,int14p,ltor,litr,int14,remltor,remitor,remint14,rematom)
    
    implicit none

    include 'mpif.h'
    
    integer :: iproc,nproc,kprint,ntap,ltors,litor,int14p
    integer :: ltor(4,ltors),litr(4,litor),int14(2,int14p)
    logical :: remltor(ltors),remitor(litor),remint14(int14p),rematom(ntap)
    integer :: i,i1,j,itmp
    integer :: nerrors,ierr,iproc0,nproc_eff
    real(8) :: dum1,dum2,dum3
    real(8) :: buffer1(0:maxrep),buffer2(0:maxrep),buffer3(0:maxrep)
    CHARACTER*27 :: err_unr(4)
    CHARACTER*15 :: err_fnf
    character*80 :: errmsg
    logical :: exist
    
    ! set error messages
    
    err_unr(1)='Unrecognized command  ---> '
    err_unr(2)='Unrecognized subcommand -> '
    err_unr(3)='Unrecognized keyword ----> '
    err_unr(4)='UNSUPPORTED  COMMAND ----> '
    err_fnf = ' file not found' 
    nerrors = 0
    
    
    ! cdist = 0: restart from an old run
    ! cdist = 1: new run with a geometric progression for temperatures 
    ! cdist = 2: new run with a custom progression 
    
    cdist_select: select case(cdist)
       
    case(0)
       
       ! ****************************************************************************************
       ! ********************************* restart an old run ***********************************
       ! ****************************************************************************************
       
       ! look for REM.set; then restart (no problems if old run was with batteries; REM set is OK)
       inquire(FILE="REM.set",EXIST=exist)
       if(exist) then
          open(unit=999,file="REM.set")
          read(999,*)
          read(999,*)
          do i = 1,nproc
             read(999,*,end=1001) itmp, dum1, dum2, dum3
             if(para_index.eq.i) then 
                rem_factor(1) = dum1
                rem_factor(2) = dum2
                rem_factor(3) = dum3
             endif
          enddo
          close(999)
          
          goto 1002
1001      errmsg = 'End of file encountered in REM.set'
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
          
1002      write(*,8099) para_index,(rem_factor(i),i=1,3)
8099      FORMAT( "-- FROM RESTART FILE:", /, "-- REPLICA INDEX",I10,/,"-- SCALING_FACTOR",3F10.4/) 
       else
          errmsg='Cannot find REM.set for restart'
          CALL xerror(errmsg,80,1,30)
          nerrors=nerrors+1
       endif
       
       ! assign to each process a replica index and a scaling 
       
    case(1)
       
       ! ****************************************************************************************
       ! *********************************  run from scratch  ***********************************
       ! ******************************  geometric progression  *********************************
       ! ****************************************************************************************
       
       para_index = iproc + 1 
       
       ! use a geometric progression
       do i = 1,3
          if(nproc.gt.nbatteries) then 
             if(mod(nproc,nbatteries).ne.0) THEN 
                if(iproc.eq.0) write(*,12754) nproc, nbatteries
12754           FORMAT( '************************************************' / &
               '*          FATAL ERROR!!!                      *' / &
               '* -nproc is not a multiple of nbatteries-      *' / &
               '* nproc   =',i4,' nbatteries =',i4,'           *' / &
               '* Action: change nbatteries in &READ_REM       *' / &
               '* or change nproc in mpiexec such that         *' / &
               '* NPROC = n*nbatteries with n integer          *' / &
               '************************************************' /)
                CALL MPI_Finalize(ierr)
                STOP
             end if
             if(nbatteries.gt.1) THEN ! produces shark teeth scaling protocol for batteries 
                nproc_eff = nproc/nbatteries 
                iproc0=iproc-int(iproc/nproc_eff)*nproc_eff
!               write(*,*) "from REM", iproc0,iproc,nproc_eff
                rem_factor(i) = rem_factor_max(i)**(dfloat(iproc0)/(nproc_eff -1.d0))
             else
                rem_factor(i) = rem_factor_max(i)**(dfloat(iproc)/(nproc -1.d0))
             end if
          else                
             ! if the REM run has only one replica, the scaling factors are the final ones
             rem_factor(i) = rem_factor_max(i)
          endif
       enddo
       
    case(2)
       
       ! ****************************************************************************************
       ! *********************************  run from scratch  ***********************************
       ! ***********************************  read a file  **********************************
       ! ****************************************************************************************
       
       para_index = iproc + 1 
       
       ! read scaling factors from a file
       INQUIRE(FILE="../REM.set",EXIST=exist)
       IF(exist) THEN
          write(kprint,*) 'Found a REM.set file'
          OPEN(unit=999,file='../REM.set')
          read(999,*)
          read(999,*)
          do i = 0,nproc-1
             read(999,*,end=1003) itmp, dum1, dum2, dum3
             if(i.eq.iproc) then 
                rem_factor(1) = dum1
                rem_factor(2) = dum2
                rem_factor(3) = dum3
             endif
          enddo
          CLOSE(999)
          
          goto 1004
1003      errmsg = 'End of file encountered in REM.set'
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
1004      continue
          
       ELSE
          nerrors = nerrors + 1
          errmsg = 'REM.set' // err_fnf 
          CALL xerror(errmsg,80,1,30)
       ENDIF
       
    case(3)
       
       ! ****************************************************************************************
       ! *********************************  run from old pdb  ***********************************
       ! ***********************************               **********************************
       ! ****************************************************************************************
       
       ! read scaling factors from a file
       INQUIRE(FILE="../REM.set",EXIST=exist)
       IF(exist) THEN
          write(kprint,*) 'Found a REM.set file'
          OPEN(unit=999,file='../REM.set')
          read(999,*)
          read(999,*)
          do i = 1,nproc
             read(999,*,end=1001) itmp, dum1, dum2, dum3
             if(para_index.eq.i) then 
                rem_factor(1) = dum1
                rem_factor(2) = dum2
                rem_factor(3) = dum3
             endif
          enddo
          CLOSE(999)
          
          goto 2004
2003      errmsg = 'End of file encountered in REM.set'
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
2004      continue
          
       ELSE
          nerrors = nerrors + 1
          errmsg = 'REM.set' // err_fnf 
          CALL xerror(errmsg,80,1,30)
       ENDIF

    case default
       
       nerrors = nerrors + 1
       errmsg = err_unr(3) // 'in REM SETUP: must be 0, 1 or 2'
       CALL xerror(errmsg,80,1,30)
       
    end select cdist_select
    
    ! if errors were found, stop
    
    if(nerrors.gt.0) then 
       errmsg= ' ERRORS WHILE VERIFYING REM INPUT' 
       call xerror(errmsg,80,1,2)
       stop
    else
       WRITE(kprint,'(5x,a)') 'REM INPUT OK!!'
    endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 
    
    ! ****************************************************************************************
    ! ****************  if all is ok, exchange informations among replicas *******************
    ! **************************************************************************************** 
    
    ! build a table for the replica indexes
    call MPI_ALLGATHER(para_index,1,MPI_INTEGER,rem_table,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    
    ! the root process gather the scaling factors from all the replicas
    call MPI_ALLGATHER(rem_factor(1),1,MPI_DOUBLE_PRECISION,buffer1 ,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    
    ! the root process gather the scaling factors from all the replicas
    call MPI_ALLGATHER(rem_factor(2),1,MPI_DOUBLE_PRECISION,buffer2,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    
    ! the root process gather the scaling factors from all the replicas
    call MPI_ALLGATHER(rem_factor(3),1,MPI_DOUBLE_PRECISION,buffer3,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    
    ! build rem_mat
    do i = 0,nproc-1
       rem_mat(rem_table(i),1) = buffer1(i)
       rem_mat(rem_table(i),2) = buffer2(i)
       rem_mat(rem_table(i),3) = buffer3(i)
    enddo
    
    if(cdist .GT. 0) then 
       OPEN(unit=999,file="REM.set")
       write(999,'(a)') 'Scaling factors:'
       write(999,'(a)') &
            'n.Ensemble      Bending+Bonding               Torsions+1-4'// &
            '                  Non-Bonded'
       do i = 1, nproc
          write(999,'(i6,1x,3g30.15)') i,( rem_mat(i,j), j=1, 3 )
       enddo
       CLOSE(999)
    endif
    
    ! ****************************************************************************************
    ! **************************  prepare a solute tempering run  ****************************
    ! ****************************************************************************************
    
    if(rem_segment) then 
       
       ! build three logical arrays, remltor and remint14 and rematom
       remltor  = .FALSE.
       remitor  = .FALSE.
       remint14 = .FALSE.
       rematom  = .FALSE.
       
       do i = 1,ltors
          do j = 1,rem_groups
             do i1 = 1,4
                ! if the atoms of the torsion overlap with the subsets of atoms, then 
                if((ltor(i1,i).GE.rem_group(1,j)).AND.(ltor(i1,i).LE.rem_group(2,j))) remltor(i) = .TRUE.
             enddo
          enddo
          if(remltor(i)) then 
             ! the 14 interaction will be scaled too
             ! find the corresponding 14 interaction
             do j = 1,int14p
                if((int14(1,j).EQ.ltor(1,i)).AND.(int14(2,j).EQ.ltor(4,i))) then
                   remint14(j) = .TRUE.
                endif
             enddo
          endif
       enddo

       do i = 1,litor
          do j = 1,rem_groups
             do i1 = 1,4
                ! if the atoms of the torsion overlap with the subsets of atoms, then 
                if((litr(i1,i).GE.rem_group(1,j)).AND.(litr(i1,i).LE.rem_group(2,j))) remitor(i) = .TRUE.
             enddo
          enddo
       enddo
       

       do i = 1,ntap
          do j = 1,rem_groups
             if((i.GE.rem_group(1,j)).AND.(i.LE.rem_group(2,j))) rematom(i) = .TRUE.
          enddo
       enddo

    endif
    
    step_rem = 0
    
  end subroutine init_rem

!-----------------------------------------------------------

  
  subroutine rem_exchange(time,nrject,nstep,iproc,nproc,t)
    !-------------------------------------------------------
    !
    ! performs replica exchanges
    !
    !-------------------------------------------------------
    !
    ! para_index: the order in temperature
    ! rem_factor: the scaling factor value
    ! rem_pot: its potential
    !
    !-------------------------------------------------------
    ! Four replica example
    !
    ! Even number of REM steps           ! Odd number of REM steps
    ! ------------------- para_index = 4 !  ------------------- para_index = 4
    !         | |                        !
    ! ------------------- para_index = 3 !  ------------------- para_index = 3
    !                                    !         | |
    ! ------------------- para_index = 2 !  ------------------- para_index = 2
    !         | |                        !
    ! ------------------- para_index = 1 !  ------------------- para_index = 1
    !
    !-------------------------------------------------------
    
    use unit
    implicit none

    include 'mpif.h'

    integer, intent(in) :: nstep
    integer, intent(in) :: nrject
    integer, intent(in) :: iproc,nproc
    real(8), intent(in) :: time
    real(8), intent(in) :: t
    ! local variables 
    integer :: i,j
    integer :: ierr
    integer :: low,high
    integer :: para_index1
    integer :: status(MPI_STATUS_SIZE)
    real(8) :: rem_delta,rem_boltz
    real(8) :: random_number
    real(8) :: ranf
    real(8) :: rem_pot1(potterms)
    external ranf

    ! attempts exchanges only after thermalization is over.
    if(nstep.gt.nrject) THEN  
       if(mod(time*nstep,rem_ts).lt.1.d-10) then 
          step_rem = step_rem + 1
          ! initialize low and high
          low = 999             
          high  = 999   
          ! the root process generates a random number random_number
          if(iproc.eq.0) random_number = ranf(0.d0) 
          ! MPI_BCAST broadcasts the random number from the root process to all the processes
          call MPI_BCAST(random_number,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
          if(mod(step_rem+para_index,2).eq.0) then
             ! Find the partner replica and assign the roles in the pair: 
             ! the low replica will take all the decisions  
             if(para_index.ne.1) then  
                high = iproc
                do i = 0,nproc-1
                   if(rem_table(i).eq.(para_index-1)) low = i
                enddo
             endif
          else
             if(para_index.ne.nproc) then
                low = iproc
                do i = 0,nproc-1
                   if(rem_table(i).eq.(para_index+1)) high = i
                enddo
             endif
          endif
          
          ! The high replica send its potential to the low replica
          if(iproc.eq.high) then 
             call MPI_SEND(rem_pot,3,MPI_DOUBLE_PRECISION,low,0,MPI_COMM_WORLD,ierr)
          elseif(iproc.eq.low) then 
             call MPI_RECV(rem_pot1,3,MPI_DOUBLE_PRECISION,high,0,MPI_COMM_WORLD,status,ierr)
          endif
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)            
          
          if(iproc.eq.low) then 
             rem_delta = 0.d0
             do i = 1,3
                rem_delta = rem_delta + &
                     (rem_mat(para_index+1,i)-rem_mat(para_index,i)) * &
                     (rem_pot1(i)-rem_pot(i))
             enddo
             rem_delta = rem_delta / (gascon*t*0.001d0)
             if(rem_delta.ge.0.d0) then 
                rep_acc(para_index) = rep_acc(para_index) + 1
                para_index1 = para_index
                para_index = para_index + 1
                do i = 1,3
                   rem_factor(i) = rem_mat(para_index,i)
                enddo
             else
                rem_boltz = dexp(rem_delta)
                if(random_number.le.rem_boltz) then ! accept the exchange
                   rep_acc(para_index) = rep_acc(para_index) + 1
                   para_index1 = para_index
                   para_index = para_index + 1
                   do i = 1,3
                      rem_factor(i) = rem_mat(para_index,i)
                   enddo
                else              
                   ! reject the exchange
                   ! the index of the high replica remains the same
                   para_index1 = para_index + 1 
                endif
             endif
             ! the low replica sends the new replica index to the high replica
             call MPI_SEND(para_index1,1,MPI_INTEGER,high,1 ,MPI_COMM_WORLD,ierr)
          elseif(iproc.eq.high) then 
             call MPI_RECV(para_index1,1,MPI_INTEGER,low,1,MPI_COMM_WORLD,status,ierr)
             para_index = para_index1
             do i = 1,3
                rem_factor(i) = rem_mat(para_index,i)
             enddo
          endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)            
          
          ! update all the vectors
          call MPI_ALLGATHER(para_index,1,MPI_INTEGER,rem_table,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
          
          if(mod(time*nstep,rem_print).lt.1.d-10.and.nstep.gt.1) then 
             call MPI_ALLREDUCE (rep_acc,rem_acc, nproc,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,ierr)
             if(iproc.eq.0) then
                do i = 1,nproc-1
                   write(*,'(15x,i5,a7,i5,a9,f8.3,5x,i10)') &
                        i, ' < = > ',  i+1,' Nacc/N% ', 2.*float(rem_acc(i))*100./float(step_rem),rem_acc(i)
                enddo
             endif
             if(mod(time*nstep,rem_print*2).lt.1.d-10.and.time*nstep.lt.rem_adjust_ratio) THEN
                if(iproc.eq.0) write(6,*) " Balancing acceptance ratio.... "
                call rem_balance(nproc,iproc)
                if(iproc.eq.0) write(6,*) " Rewriting the REM.set .." 
             endif
             !all processes rewrite updated REM.set 
             if(mod(time*nstep,rem_print*2).lt.1.d-10.and.time*nstep.lt.rem_adjust_ratio) THEN
                step_rem=0
                do i=1,nproc
                   rep_acc(i)=0
                end do
                do i=1,nproc-1
                   rem_acc(i)=0
                end do
                CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
                OPEN(unit=999,file="REM.set")
                write(999,'(a)') 'Scaling factors:'
                write(999,'(a)') &
                     'n.Ensemble      Bending+Bonding               Torsions+1-4'// &
                     '                  Non-Bonded'
                do i = 1, nproc
                   write(999,'(i6,1x,3g30.15)') i,( rem_mat(i,j), j=1, 3 )
                enddo
                CLOSE(999)
             endif
          endif
       endif 
    endif

  end subroutine rem_exchange

  SUBROUTINE rem_balance(nproc,iproc)

    implicit none

    include 'mpif.h'

    integer :: nproc,iproc
    ! global (module) variables
    ! rem_mat (nrep,3) REM.set matrix

    ! local variables 
    integer :: i,j,k,n1
    integer :: ierr
    integer :: nproc_eff
    real*8 :: ratio(nproc),dn(0:nproc,3),d(0:nproc,3), mean_ratio,sig,c,old(3),rnew(3)

    ! computes average exchange ratio across batteries
    nproc_eff=nproc/nbatteries
    k=0
    c=2.d0   ! initial (and hopefully final) value of c
    d=0
    dn=0
    ratio=0.d0
    do i=1,nproc-1
       k=k+1
       ratio(k)=2.*float(rem_acc(i))/float(step_rem)+ratio(k)
       if(k.eq.nproc_eff) k=0
    end do
    ratio=ratio/float(nbatteries)
    mean_ratio=0.d0
    sig=0.0d0
    n1=nproc_eff-1
    do i=1,n1
       mean_ratio=mean_ratio+ratio(i)
       sig=sig+ratio(i)**2
    end do
    mean_ratio=mean_ratio/dfloat(n1)
    sig=dsqrt(sig/dfloat(n1)-mean_ratio**2)
    if(iproc.eq.0) write(6,'(A9,f5.2,A10,f6.3)') "  Mean = ",mean_ratio, "   Sigma =",sig
    ! compute difference in the REM set (torsional potential as referrence)
    old=1.0
    do j=1,3
       do i=2,nproc_eff
          d(i-1,j)=old(j)-rem_mat(i,j) 
          old(j)=rem_mat(i,j)
       end do
    end do
    ! update differences using balancing formula based on shift from mean acc ratio
100 continue
    do j=1,3
       do i=1,n1
          dn(i,j)=0.d0
          ! change only active scaling factors  
          if(abs(d(i,j)).gt.1D-15) THEN
             dn(i,j)=d(i,j)+c*sig*(ratio(i)-mean_ratio)
          END IF
       end do
    end do
    
!    do i=1,n1
!       if(iproc.eq.0) write(6,"(i5,2(3f10.3,3x))") i,(dn(i,j),j=1,3),(d(i,j),j=1,3)
!    end do

    ! update basic rem_mat using torsional potential 
  
    rnew=1.d0
    do i=0,nproc_eff-1
       do j=1,3
          rem_mat(i+1,j)=rnew(j)-dn(i,j)
          if(rem_mat(i+1,j).lt.0) THEN
             c=1.d0
             if(iproc.eq.0) write(6,*) "Negative scaling: retrying balance with c=1 ..."
             goto 100
          ENDIF
          rnew(j)=rnew(j)-dn(i,j)
       end do
    end do
    ! propagate to all batteries and rewrite the REM.set file
    k=0
    if(nbatteries.gt.1) THEN  
       do i=nproc_eff+1,nproc
          k=k+1
          do j=1,3
             rem_mat(i,j)=rem_mat(k,j)
          end do
          if(k.eq.nproc_eff) k=0
       end do
    ENDIF
  END SUBROUTINE rem_balance

#endif

  LOGICAL FUNCTION seg_select(kind,rematom,a,b)
    IMPLICIT NONE
    INTEGER a,b
    INTEGER kind
    LOGICAL rematom(*)
    
    if(kind.EQ.0) then 
       seg_select = rematom(a).or.rematom(b)
    elseif(kind.EQ.1) then 
       seg_select = rematom(a).and.rematom(b)
    elseif(kind.EQ.2) then 
       seg_select = rematom(a).neqv.rematom(b)
    endif
    
  END FUNCTION seg_select

  
end module rem
