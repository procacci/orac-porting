!********************************************************************
!                                                                   *
!   Reference author for this part of code: Riccardo Chelli         *
!   Department of Chemistry, University of Firenze, Italy           *
!   E-mail: riccardo.chelli@unifi.it                                *
!   Module introduced for Serial Generalized Ensemble simulations   *
!   Last update: 6 June 2012                                        *
!                                                                   *
!********************************************************************

module sge

  implicit none
  save

  integer :: sge_df_opt        ! 0: use fixed free energies in SGE simulation; 1:update free energies in SGE simulation (using BAR-SGE)
  integer :: nstates           ! number of generalized-ensemble states in SGE simulation
  integer :: sge_update        ! number of Df estimates used to update free energies in SGE simulation
  integer :: sge_update_rst    ! number of Df estimates used to update free energies in SGE simulation (in restart file)
  integer :: sge_update_fly    ! number of Df estimates used to update free energies in SGE simulation (for change on the fly)
  integer :: nstep_fly         ! tag for changing Df on the fly in SGE simulation
  integer :: nfly_df           ! counter for Df change during the SGE simulation
  integer :: ncycle_reject     ! counter for cyclic SGE  (rej cycles)
  integer :: ncycle            ! counter  for cyclic SGE  (all cycles) 
  integer :: target_replica    ! target replica in cyclic SGE
  integer :: nfreqsge          ! printout integer (active when always_accept=T)
  integer :: print_always      ! printout integer (1:pdb; 2:xyz; active when always_accept=T)
  INTEGER :: write_sge_work    ! print work (always_accept); if (d) print it at all times.
  integer :: ud                ! variable used in the upward/downward replica transition scheme
  real(8) :: sge_tsave_work    ! work dumping interval in SGE simulation
  real(8) :: sge_tsave_pmf     ! potential of mean force updating interval in SGE simulation
  real(8) :: sge_ts            ! exchange time interval in SGE simulation
  real(8) :: sge_print         ! time interval for printing acceptance ratio
  real(8) :: sge_f             ! saving time interval for SGE_WHAM file (replica, unscaled energies/potential parameters)
  real(8) :: wtot,time_fs
  logical :: sim_tempering     ! if true, a SGE simulation is performed
  logical :: sim_tempering_sge ! if true, a SGE simulation in the collective coordinate space is performed
  logical :: zero_free_energy  ! if true, the free energies are zeroed when a run is restarted
  character(3) :: sge_transition_scheme ! scheme for replica transitions (SEO or DEO)

  integer, allocatable :: estimate_df(:), estimate_df_up(:), estimate_df_down(:)
  integer, allocatable :: sge_up(:), sge_up_jar(:), sge_down(:), sge_down_jar(:)
  integer, allocatable :: nconst_down(:), nw_down_tot_old(:)
  integer, allocatable :: nconst_up(:), nw_up_tot_old(:), cnt_upd(:)
  real(8), allocatable :: w_up(:,:), w_up_jar(:,:), w_down(:,:), w_down_jar(:,:)
  real(8), allocatable :: df(:), dfi(:), df_var(:), df_av(:), ndf_av(:), ndf_av_fly(:), &
       df_down(:), df_up(:), df_fixed(:), dfi_inst(:,:), df_var_inst(:,:)
  real(8), allocatable :: acc_up(:), acc_down(:), tot_up(:), tot_down(:)
  real(8), allocatable :: histog(:), histog_t(:)
  real(8), allocatable :: sge_eqdist(:,:), sge_eqang(:,:), sge_eqdied(:,:)
  integer, allocatable :: idx_bonds(:), idx_bends(:), idx_tors(:)
  integer :: nbonds_added_st, nbends_added_st, nitors_added_st
  integer :: ierror, current_dim_up, current_dim_up_jar
  integer :: current_dim_down, current_dim_down_jar
  integer :: nt_points_w, ksge
  logical :: sge_bonds, sge_bends, sge_tors, lfirst_sge, sge_print_wham,go_save,use_saved,always_accept
  character(80) :: df_file

! Auxiliary variables (employed into readrs subroutine)
  integer :: nstates_rs

! Real kinds
  integer, parameter :: kr4 = selected_real_kind(6,37)     ! Single precision real 
  integer, parameter :: kr8 = selected_real_kind(15,307)   ! Double precision real 
  integer, parameter :: kr16 = selected_real_kind(30,1000) ! Quadruple precision real 
! Integer kinds 
  integer, parameter :: ki4 = selected_int_kind(9)  ! Single precision integer 
  integer, parameter :: ki8 = selected_int_kind(18) ! Double precision integer 
! Complex kinds 
  integer, parameter :: kc4 = kr4 ! Single precision complex 
  integer, parameter :: kc8 = kr8 ! Double precision complex 


contains

!===============  SGE_DRIVE  =======================================================
  subroutine sge_drive ( bond_added, bend_added, tors_added, radfact, nstep, iproc, &
     ! Energies
       uslvbon, uslvben, uslvitor, ubond, ubend, uitors, efact, uslvtor, conf_bnd_slv_n1, coul_bnd_slv_n1, uptors, &
       conf_bnd_slt_n1, coul_bnd_slt_n1, ucns_h, ucos_h, ucns_l, ucos_l, ucns_m, ucos_m, urcs_h, urcs_l, urcs_m,   &
       coul_bnd_slv, self_slv, urcp_h, urcp_l, urcp_m, ucop_h, ucop_l, ucop_m, ucnp_h, ucnp_l, ucnp_m,             &
       coul_bnd_slt, self_slt, urcsp_h, urcsp_l, urcsp_m, ucosp_h, ucosp_l, ucosp_m, ucnsp_h, ucnsp_l, ucnsp_m,    &
       eer_m, eer_l, eer_h, &
     ! Other energy related variables
       hpot, pext, pucek, uceh, ucek, ucepr, volume, cpress )

    use cpropar, only: eqdist, eqang, eqdied, lstretch, lbend, litor, nbonds_added, nbends_added, nitors_added, &
         force_cost, force_ang, force_died, potbo, potbe, potit, nrject, t, time
    use unit, only: gascon,kprint,kwork
    use rem, only: para_index, rem_pot, rem_mat, rem_segment, rem_hnb_energy, rem_lnb_energy, rem_mnb_energy, rem_n1_energy

    implicit none

    real(8), external :: ranf
    integer, intent(in) :: nstep, iproc
    real(8), intent(in) :: radfact, bond_added(*), bend_added(*), tors_added(*)
    real(8), intent(in) :: uslvbon, uslvben, uslvitor, ubond, ubend, uitors, efact, uslvtor, conf_bnd_slv_n1, coul_bnd_slv_n1, &
         uptors, conf_bnd_slt_n1, coul_bnd_slt_n1, ucns_h, ucos_h, ucns_l, ucos_l, ucns_m, ucos_m, urcs_h, urcs_l, urcs_m,     &
         coul_bnd_slv, self_slv, urcp_h, urcp_l, urcp_m, ucop_h, ucop_l, ucop_m, ucnp_h, ucnp_l, ucnp_m, coul_bnd_slt,         &
         self_slt, urcsp_h, urcsp_l, urcsp_m, ucosp_h, ucosp_l, ucosp_m, ucnsp_h, ucnsp_l, ucnsp_m, eer_m, eer_l, eer_h
    real(8), intent(in) :: hpot, pext, pucek, uceh, ucek, volume, ucepr
    logical, intent(in) :: cpress
    real(8) :: n1_energy, un0, un1, unb, rest_n1, rest_nb, ene_nb
    real(8) :: acc_down_t, acc_up_t, tot_down_t, tot_up_t
    real(8) :: dff, populat, random_number,sge_delta
    integer :: i, idx
    character(80) :: df_fly_file
    logical :: exist_df
    integer :: ios

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE), ierr
#endif

! Current simulation time (in fs units)
    time_fs = time * nstep

!---------------------------------------------------------------------------------------------------
! Initialization of the SGE_ENERGY, SGE_DF and SGE_WHAM files
    if ( lfirst_sge ) THEN 
       if(always_accept) THEN 
          wtot=0.d0
          ncycle=0 
          ncycle_reject=0 
       END IF
       call init_files ( iproc, cpress )
    end if

!---------------------------------------------------------------------------------------------------
! Update rem_pot if necessary and calculate the unscaled potential energy
    if ( sim_tempering_sge ) then
       if ( mod(time_fs,sge_ts) < 1.d-10 ) then
          un0 = ( uslvbon + uslvben + uslvitor + ubond + ubend + uitors ) * efact / 1000.d0
          un1 = ( uslvtor + conf_bnd_slv_n1 + coul_bnd_slv_n1 + uptors + conf_bnd_slt_n1 + coul_bnd_slt_n1 ) * efact / 1000.
          unb = ( ucns_h + ucos_h + ucns_l + ucos_l + ucns_m + ucos_m + urcs_h + urcs_l + urcs_m + coul_bnd_slv + self_slv + &
               urcp_h + urcp_l + urcp_m + ucop_h + ucop_l + ucop_m + ucnp_h + ucnp_l + ucnp_m + coul_bnd_slt + self_slt +    &
               urcsp_h + urcsp_l + urcsp_m + ucosp_h + ucosp_l + ucosp_m + ucnsp_h + ucnsp_l + ucnsp_m + eer_m + eer_l +     &
               eer_h ) * efact / 1000.d0
       endif
    else
       if ( ( mod(time_fs,sge_ts) < 1.d-10 ) .or.                         &
            ( mod(time_fs,sge_tsave_work) < 1.d-10 .and. nstep > 0 ) .or. &
            ( mod(time_fs,sge_f) < 1.d-10 ) ) then
          un0 = ( uslvbon + uslvben + uslvitor + ubond + ubend + uitors ) * efact / 1000.d0
          un1 = ( uslvtor + conf_bnd_slv_n1 + coul_bnd_slv_n1 + uptors + conf_bnd_slt_n1 + coul_bnd_slt_n1 ) * efact / 1000.
          unb = ( ucns_h + ucos_h + ucns_l + ucos_l + ucns_m + ucos_m + urcs_h + urcs_l + urcs_m + coul_bnd_slv + self_slv + &
               urcp_h + urcp_l + urcp_m + ucop_h + ucop_l + ucop_m + ucnp_h + ucnp_l + ucnp_m + coul_bnd_slt + self_slt +    &
               urcsp_h + urcsp_l + urcsp_m + ucosp_h + ucosp_l + ucosp_m + ucnsp_h + ucnsp_l + ucnsp_m + eer_m + eer_l +     &
               eer_h ) * efact / 1000.d0

          rem_pot(1) = un0

          if ( rem_segment ) then
! Calculate the difference between the total potential and the one from the subset
             n1_energy = rem_n1_energy * efact / 1000.d0
             ene_nb = ( rem_hnb_energy + rem_lnb_energy + rem_mnb_energy ) * efact / 1000.d0
             rest_n1 = un1 - n1_energy
             rest_nb = unb - ene_nb
! Only the potential terms involving the selected atoms
             rem_pot(2) = n1_energy
             rem_pot(3) = ene_nb
          else
             rest_n1 = 0.d0
             rest_nb = 0.d0
! Use all the potential
             rem_pot(2) = un1
             rem_pot(3) = unb
          endif
       endif
    endif

!------------------------------------------------------------------------------------
! Print SGE_WHAM for WHAM application
    if(.not.always_accept) THEN 
       if ( mod(time_fs,sge_f) < 1.d-10 ) then
          call print_sge_wham ( time_fs, tors_added, eqdied, bend_added, eqang, bond_added, eqdist, radfact )
       endif
    endif

!------------------------------------------------------------------------------------
! Print SGE_ENERGY
    if ( mod(time_fs,sge_ts) < 1.d-10 ) then
       call print_sge_energy ( time_fs, ucek, pucek, ucepr, uceh, pext, volume, hpot, efact, un0, un1, unb, &
            rest_n1, rest_nb, cpress )
    endif

!------------------------------------------------------------------------------------
! In rejecting phase, the SGE algorithm is skipped
    if ( nstep <= nrject ) goto 150

!-------------------------------------------------------------------------------------
! Update the arrays containing the generalized work needed to determine the upward and downward replica transitions
    if ( mod(time_fs,sge_tsave_work) < 1.d-10 .and. nstep > 0 ) then

!-------------------------------------------------------------------------------------
       call sge_size_update ( iproc )

! UPWARD REPLICA TRANSITIONS
       if ( para_index /= nstates ) then

          sge_delta = 0.d0
! Part related to MultiWindow tempering simulation
          if ( sim_tempering_sge ) then
             do i = 1, nbonds_added_st
                idx = idx_bonds(i)
                sge_delta = sge_delta + force_cost(idx) *                    &
                     ( ( bond_added(idx) - sge_eqdist(para_index+1,i) )**2 - &
                       ( bond_added(idx) - sge_eqdist(para_index,  i) )**2 )
             enddo
             do i = 1, nbends_added_st
                idx = idx_bends(i)
                sge_delta = sge_delta + force_ang(idx) *                            &
                     ( ( bend_added(idx) - sge_eqang(para_index+1,i)*radfact )**2 - &
                       ( bend_added(idx) - sge_eqang(para_index,  i)*radfact )**2 )
             enddo
             do i = 1, nitors_added_st
                idx = idx_tors(i)
                sge_delta = sge_delta + force_died(idx) *                            &
                     ( ( tors_added(idx) - sge_eqdied(para_index+1,i)*radfact )**2 - &
                       ( tors_added(idx) - sge_eqdied(para_index,  i)*radfact )**2 )
             enddo
! Conversion from kcal/mol to kJ/mol
             sge_delta = sge_delta * 4.184d0
! Part related to SOLUTE TEMPERING
          else
             do i = 1, 3
                sge_delta = sge_delta + ( rem_mat(para_index+1,i) - rem_mat(para_index,i) ) * rem_pot(i)
             enddo
          endif

          sge_delta = sge_delta / ( gascon * t * 0.001d0 )
          sge_up(para_index) = sge_up(para_index) + 1
          w_up(sge_up(para_index),para_index) = sge_delta

          if ( estimate_df(para_index) == 0 ) then
             sge_up_jar(para_index) = sge_up_jar(para_index) + 1
             w_up_jar(sge_up_jar(para_index),para_index) = sge_delta
          endif
       endif

! DOWNWARD REPLICA TRANSITIONS
       if ( para_index /= 1 ) then

          sge_delta = 0.d0
! Part related to Multiwindow Tempering simulation
          if ( sim_tempering_sge ) then
             do i = 1, nbonds_added_st
                idx = idx_bonds(i)
                sge_delta = sge_delta + force_cost(idx) *                   &
                     ( ( bond_added(idx) - sge_eqdist(para_index-1,i))**2 - &
                       ( bond_added(idx) - sge_eqdist(para_index,  i))**2 )
             enddo
             do i = 1, nbends_added_st
                idx = idx_bends(i)
                sge_delta = sge_delta + force_ang(idx) *                            &
                     ( ( bend_added(idx) - sge_eqang(para_index-1,i)*radfact )**2 - &
                       ( bend_added(idx) - sge_eqang(para_index,  i)*radfact )**2 )
             enddo
             do i = 1, nitors_added_st
                idx = idx_tors(i)
                sge_delta = sge_delta + force_died(idx) *                            &
                     ( ( tors_added(idx) - sge_eqdied(para_index-1,i)*radfact )**2 - &
                       ( tors_added(idx) - sge_eqdied(para_index,  i)*radfact )**2 )
             enddo
! Conversion from kcal/mol to kJ/mol
             sge_delta = sge_delta * 4.184d0
! Part related to SOLUTE TEMPERING
          else
             do i = 1, 3
                sge_delta = sge_delta + ( rem_mat(para_index-1,i) - rem_mat(para_index,i) ) * rem_pot(i)
             enddo
          endif

          sge_delta = sge_delta / ( gascon * t * 0.001d0 )
          sge_down(para_index) = sge_down(para_index) + 1
          w_down(sge_down(para_index),para_index) = sge_delta

          if ( estimate_df(para_index-1) == 0 ) then
             sge_down_jar(para_index) = sge_down_jar(para_index) + 1
             w_down_jar(sge_down_jar(para_index),para_index) = sge_delta
          endif
       endif

    endif

!-----------------------------------------------------------------------------
! 1) Update of free energy differences
! 2) Print (dimensionless) free energy difference (DF) and related error 
! 3) Print histogram of ensemble-population for diagnostic
    if ( ( mod(time_fs,sge_tsave_pmf) < 1.d-10 .and. nstep > 0 ) .or. nstep_fly == 0 ) then

#ifdef _MPI_
       df_fly_file = '../SGE_DF_FLY.set'
#else
       df_fly_file = 'SGE_DF_FLY.set'
#endif
       inquire ( file=df_fly_file, exist=exist_df )
       sge_update_fly = 0

       if ( exist_df ) then
          open( unit=666, file=df_fly_file, form='formatted' )
          read(666,*,end=700) sge_update_fly
700       close(666)
          if ( sge_update_fly > 0 ) then
             if ( iproc == 0 ) then
#ifdef _MPI_
                df_fly_file = '../SGE_DF_FLY.dat'
#else
                df_fly_file = 'SGE_DF_FLY.dat'
#endif
                if ( nstep_fly == 0 ) open( unit=444, file=df_fly_file, form='formatted' )
                rewind(444)

                nfly_df = 0
                do
                   read(444,*,iostat=ios)
                   write(*,*) "main: nfly_df, ios= ", nfly_df, ios
                   if(ios<0) exit
                   nfly_df = nfly_df + 1
                enddo

#if  ( __GNUC__ < 4 ||  (__GNUC__ == 4 && __GNUC_MINOR__ < 6 )) || defined AIX
                ! --- backspace of file 444 skipped on gfor<4.6 and AIX ---
#else 
                backspace(444)
#endif

                if ( nstep_fly == 0 ) then
                   write(*,'(a,1x,i10)') 'Initial estimate of weighted '//  &
                        'free energies (available number of samples from file SGE_DF_FLY.dat):', nfly_df
                   if ( nfly_df > 0 ) then
                      do i = 1, nstates - 1
                         call sge_df_update_2 ( i )
                      enddo
                   endif
                endif
             endif
#ifdef _MPI_
             if ( nstep_fly == 0 ) then
                do i = 1, nstates - 1
                   call mpi_bcast( df(i), 1, mpi_double_precision, 0, mpi_comm_world, ierr )
                   call mpi_bcast( ndf_av_fly(i), 1, mpi_double_precision, 0, mpi_comm_world, ierr )
                enddo
             endif
#endif
          endif
       endif

       if ( nstep_fly == 1 ) then
! Free energy update
          call sge_df_update ( time_fs, iproc )

! Writing DF, average DF and error
          if ( iproc == 0 .and. (.not.always_accept)) then
             if ( sge_update_fly > 0 ) then
                write(889,'(f14.4,1000(g19.10,1x))') &
                     time_fs, ( df(i), i=1, nstates-1 ), ( sqrt(1.d0/ndf_av_fly(i)), i=1, nstates-1 )
             else
                write(889,'(f14.4,1000(g19.10,1x))') &
                     time_fs, ( df(i), i=1, nstates-1 ), ( sqrt(1.d0/ndf_av(i)), i=1, nstates-1 )
             endif
          endif

! Writing histogram of ensemble-population
          populat = 0.d0
          do i = 1, nstates
#ifdef _MPI_
             call mpi_reduce( histog(i), histog_t(i), 1, mpi_double_precision,mpi_sum, 0, mpi_comm_world, ierr )
#else
             histog_t(i) = histog(i)
#endif
             populat = populat + histog_t(i)
          enddo

          if ( iproc == 0 ) then
             open( unit=890, file='SGE_HISTOG', form='formatted' )
             write(890,'(a)') '# Ensemble   Population(%)   Overall population'
             do i = 1, nstates
                write(890,'(3x,i4,5x,f12.3,5x,g20.12)') i, histog_t(i) / populat * 100.d0, histog_t(i)
             enddo
             close(unit=890)
          endif
       endif
       nstep_fly = 1
    endif

!------------------------------------------------------------------------
! Attempiting upward and downward replica transitions
! (also update the histogram of the ensemble-population)
    if ( mod(time_fs,sge_ts) < 1.d-10 .and. nstep > 0 ) then

! Update the histogram of the ensemble-population
       histog(para_index) = histog(para_index) + 1.d0

! Applies the chosen replica transition scheme
! (SEO or DEO), see Denschlag et al. JCTC 5,2847 (2009)
       if ( sge_transition_scheme == 'SEO' ) then
! "Stochastic even/odd" (SEO) scheme
          random_number = ranf()
          if ( random_number > 0.5 ) then
             ud = 1
          elseif ( random_number < 0.5 ) then
             ud = -1
          else
             ud = 0
          endif
       else
! "Deterministic even/odd" (DEO) scheme
         if ( mod(int(time_fs/sge_ts)+para_index,2) == 0 ) then
             ud = 1
          elseif ( mod(int(time_fs/sge_ts)+para_index,2) /= 0 ) then
             ud = -1
          else
             ud = 0
          endif
       endif
       
! Upward transitions
       if ( ud == 1 .and. para_index /= nstates ) then

          if ( sge_df_opt == 1 ) then
             if ( estimate_df(para_index) == 1 ) then
                dff = df(para_index)      ! The weighted average of DF is used for the acceptance ratio
             elseif ( estimate_df_up(para_index) == 1 ) then
                dff = df_up(para_index)   ! The DF value obtained using the Jarzynski equality is used for the acceptance ratio
                                          ! (using work samples  from  para_index  to  para_index+1)
             elseif ( estimate_df_down(para_index) == 1 ) then
                dff = df_down(para_index) ! The DF value obtained using the Jarzynski equality is used for the acceptance ratio
                                          ! (using work samples  from  para_index+1  to  para_index)
             else
                dff = 0.
             endif
          else
             dff = df_fixed(para_index)
          endif

          call temp_up_down ( bond_added, bend_added, tors_added, force_cost, force_ang, force_died, eqdist,   &
               eqang, eqdied, potbo(1,2), potbe(1,2), potit(1,2), dff, 1, t, radfact, lstretch,                &
               lbend, litor, nbonds_added, nbends_added, nitors_added, iproc )

! Downward transitions
       elseif ( ud == -1 .and. para_index /= 1 ) then

          if ( sge_df_opt == 1 ) then
             if ( estimate_df(para_index-1) == 1 ) then
                dff = df(para_index-1)      ! The weighted average of DF is used for the acceptance ratio
             elseif ( estimate_df_down(para_index-1) == 1 ) then
                dff = df_down(para_index-1) ! The DF value obtained using the Jarzynski equality is used for the acceptance ratio
                                            ! (using work samples  from  para_index  to  para_index-1)
             elseif ( estimate_df_up(para_index-1) == 1 ) then
                dff = df_up(para_index-1)   ! The DF value obtained using the Jarzynski equality is used for the acceptance ratio
                                            ! (using work samples  from  para_index-1  to  para_index)
             else
                dff = 0.
             endif
          else
             dff = df_fixed(para_index-1)
          endif

          call temp_up_down ( bond_added, bend_added, tors_added, force_cost, force_ang, force_died, eqdist, eqang,    &
               eqdied, potbo(1,2), potbe(1,2), potit(1,2), dff, -1, t, radfact, lstretch, lbend, litor,                &
               nbonds_added, nbends_added, nitors_added, iproc )

       endif
    endif

!------------------------------------------------------------------------------------------------
! Print acceptance ratio for replica transitions
    if ( mod(time_fs,sge_print) < 1.d-10 .and. nstep > 0 ) then
       if(.not.always_accept) THEN 
          do i = 1, nstates-1
#ifdef _MPI_
             call mpi_reduce( acc_up(i), acc_up_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
             call mpi_reduce( acc_down(i+1), acc_down_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
             call mpi_reduce( tot_up(i), tot_up_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
             call mpi_reduce( tot_down(i+1), tot_down_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
#else
             acc_up_t = acc_up(i)
             acc_down_t = acc_down(i+1)
             tot_up_t = tot_up(i)
             tot_down_t = tot_down(i+1)
#endif
             if ( tot_up_t == 0. ) tot_up_t = 1.
             if ( tot_down_t == 0. ) tot_down_t = 1.
             if ( iproc == 0 ) then
                write(*,300) &
                     i, '-->', i+1, 'N_acc/N(%)', 100.*acc_up_t/tot_up_t, i+1, '-->', i, 'N_acc/N(%)', 100.*acc_down_t/tot_down_t
300             format(2(5x,i3,1x,a3,i3,2x,a10,f9.3))
             endif
          enddo
       ELSE 
          write(*,301) iproc,100.d0*dfloat(ncycle-ncycle_reject)/dfloat(ncycle), ncycle,ncycle_reject
301       format(" Walker: ",i8," Acceptance Ratio:",f8.2, 2i10)
       end if
    endif

150 continue

    return
  end subroutine sge_drive


!==================  READING INPUT PARAMETERS FOR SGE  =====================
  subroutine read_sge ( err_args, err_unr, err_end, err_fnf, err_open )

    use unit, only: knlist, kprint,kwork
    use rem, only: rem_segment, rem_segkind, rem_groups, cdist, rem_group, rem_factor_max

    implicit none

    integer :: iret
    character(37) :: err_args(3)
    character(20) :: err_end 
    character(27) :: err_unr(4)
    character(22) :: err_open

#ifdef _MPI_
    include 'mpif.h'
    integer :: ierr
#endif

    integer :: nword, nsevere, nwarning, j, i
    character(80) :: errmsg
    character(80) :: line, strngs(40)
    character(8) :: fmt
    character(15) :: err_fnf
    character(1) :: sep(2), comm(2)
    logical :: exist
    data sep/' ',','/comm/'(',')'/
    real(8) :: aux, dum, dum1, dum2, dum3


    rem_segment = .false.      
    sim_tempering = .true.
    zero_free_energy = .false.
    rem_groups = 0
    rem_segkind = 0
    sge_df_opt = 1

!-------------------------------------------
!  Environment parser starts here 
!-------------------------------------------

    j = 0
    nsevere = 0 
    nwarning = 0 
    line(79:80) = '  '
100 read(knlist,'(a78)',end=600) line(1:78)
    call wrenc(kprint,line)
    if (line(1:1) == '#') goto 100 
    call parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
    if (iret == 1) then 
       errmsg = 'while parsing line: toomany strings'
       call xerror(errmsg,80,1,2)
       nsevere = nsevere + 1
       goto 100
    endif

!---- Command  SETUP ----------------------------------------------
! Read the extremum value of the scaling factor
    if ( strngs(1) == 'SETUP' ) then

       if ( nword == 6 ) then
          call fndfmt(1,strngs(2),fmt)
          read(strngs(2),fmt,err=20) nstates
          call fndfmt(2,strngs(3),fmt)
          read(strngs(3),fmt,err=20) rem_factor_max(1)
          call fndfmt(2,strngs(4),fmt)
          read(strngs(4),fmt,err=20) rem_factor_max(2)
          call fndfmt(2,strngs(5),fmt)
          read(strngs(5),fmt,err=20) rem_factor_max(3)
          call fndfmt(1,strngs(6),fmt)
          read(strngs(6),fmt,err=20) cdist
       elseif ( nword == 3 ) then
          call fndfmt(1,strngs(2),fmt)
          read(strngs(2),fmt,err=20) nstates
          call fndfmt(1,strngs(3),fmt)
          read(strngs(3),fmt,err=20) cdist
          sim_tempering_sge = .true.
       else
          nsevere = nsevere + 1
          errmsg = err_args(3) //'2 or 5'
          call xerror(errmsg,80,1,30)
       endif


!------ Command  FIX_FREE_ENERGY -----------------------------------
    elseif ( strngs(1) == 'FIX_FREE_ENERGY' ) then
       sge_df_opt = 0

       if (nword == 3) then
          if ( strngs(2) == 'OPEN' ) then
             call uscrpl(strngs(3),80)
             inquire(file=strngs(3),exist=exist)
             if (exist) then
                df_file=strngs(3)
             else
                errmsg = 'Auxiliary file does not exist. Abort.'
                CALL xerror(errmsg,80,1,30)
                nsevere = nsevere + 1
             endif
          else
             errmsg = err_open
             call xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
          endif
       else
          nsevere = nsevere + 1
          errmsg = err_args(3) //'2'
          call xerror(errmsg,80,1,30)
       endif

!---- Command  ZERO_FREE_ENERGY -----------------------------
    elseif ( strngs(1) == 'ZERO_FREE_ENERGY' ) then
       zero_free_energy = .true.

!---- Command  STEP -----------------------------------------

! Read the SGE time step for replica transitions
    elseif ( strngs(1) == 'STEP' ) then
       if ( nword == 4 .or. nword == 5 ) then
          call fndfmt(2,strngs(2),fmt)
          read(strngs(2),fmt,err=20) sge_ts
          call fndfmt(2,strngs(3),fmt)
          read(strngs(3),fmt,err=20) sge_tsave_work
          call fndfmt(2,strngs(4),fmt)
          read(strngs(4),fmt,err=20) sge_tsave_pmf
          if ( nword == 5 ) then
             call fndfmt(1,strngs(5),fmt)
             read(strngs(5),fmt,err=20) sge_update
          endif
       else
          nsevere = nsevere + 1
          errmsg = err_args(3) //'3 or 4'
          call xerror(errmsg,80,1,30)
       endif

!---- Command  TRANSITION_SCHEME ---------------------------

! Read the scheme for replica transitions
    elseif ( strngs(1) == 'TRANSITION_SCHEME' ) then
       if ( nword == 2 ) then
          if ( strngs(2) /= 'SEO' .and. strngs(2) /= 'DEO' ) then
             nsevere = nsevere + 1
             errmsg = err_unr(3) // 'in TRANSITION_SCHEME(&SGE): must be "SEO" or "DEO"'
             call xerror(errmsg,80,1,30)
          endif
          sge_transition_scheme = strngs(2)
       else
          nsevere = nsevere + 1
          errmsg = err_args(3) //'1'
          call xerror(errmsg,80,1,30)
       endif
       if(always_accept) sge_transition_scheme = 'DEO'

!---- Command  SEGMENT ----------------------------------------

! Read the first and the last atom of the segment
    elseif ( strngs(1) == 'SEGMENT' ) then
       rem_segment = .true.
!------- read the line
700    read(knlist,'(a78)',end=600) line(1:78)
       call wrenc(kprint,line)
       if (line(1:1) == '#') goto 700
       call parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
       if ( strngs(1) == 'define') then
          if ( nword == 3 ) then 
             rem_groups = rem_groups + 1
             call fndfmt(1,strngs(2),fmt)
             read(strngs(2),fmt,err=20) rem_group(1,rem_groups)
             call fndfmt(1,strngs(3),fmt)
             read(strngs(3),fmt,err=20) rem_group(2,rem_groups)
          else
             nsevere = nsevere + 1
             errmsg = err_args(3) //'2'
             call xerror(errmsg,80,1,30)
          endif
       elseif ( strngs(1) == 'kind') then
          if ( nword == 2 ) then 
             if ( strngs(2) == 'intra') then
                rem_segkind = 1
             elseif ( strngs(2) == 'inter') then
                rem_segkind = 2
             endif
! Do nothing
          else
             nsevere = nsevere + 1
             errmsg = err_args(3) //'1'
             call xerror(errmsg,80,1,30)
          endif
          
       elseif ( strngs(1) == 'defaults' ) then
          continue
       elseif ( strngs(1) == ' ' ) then
          continue
       elseif ( strngs(1) == 'END' ) then
          goto 100
       else
!--- Could not find SUBCOMMAND of END
          errmsg = err_unr(2) // strngs(1)// ' or missing END'
          call xerror(errmsg,80,1,30)
          nsevere = nsevere + 1 
       endif
       goto 700

!---- Command  PRINT_ACCEPTANCE_RATIO -------------------------------------------

! read the SGE print interval for acceptance ratio
    elseif ( strngs(1) == 'PRINT_ACCEPTANCE_RATIO' ) then
       if ( nword == 2 ) then
          call fndfmt(2,strngs(2),fmt)
          read(strngs(2),fmt,err=20) sge_print
       else
          nsevere = nsevere + 1
          errmsg = err_args(3) //'1'
          call xerror(errmsg,80,1,30)
       endif

!---- Command  PRINT_WHAM ---------------------------------------

    elseif ( strngs(1) == 'PRINT_WHAM' ) then
       if ( nword == 2 ) then
          sge_print_wham = .true.
          call fndfmt(2,strngs(2),fmt)
          read(strngs(2),fmt,err=20) sge_f
          strngs(3) = 'SGE_WHAM'
          call uscrpl(strngs(3),80)
          inquire(file=strngs(3),exist=exist)
          if ( exist ) then
             call openf(ksge,strngs(3),'FORMATTED','OLD',0)
          else
             call openf(ksge,strngs(3),'FORMATTED','NEW',0)
          endif
       else
          errmsg = err_args(1)//'1'
          call xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
       endif

!---- Command  ALWAYS_ACCEPT -------------------------------------------

      ELSEIF(strngs(1).EQ. 'ALWAYS_ACCEPT' ) THEN
        always_accept=.true.
        sge_transition_scheme='DEO' 
        if(nfreqsge.eq.0) nfreqsge=1 
        if(print_always.eq.0) print_always=1
        IF(nword.eq.2) THEN 
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) target_replica
        ELSE if(nword.eq.1) THEN 
          target_replica=1
          write_sge_work=1
        ELSE if(nword.eq.3) THEN
          target_replica=1 
          write_sge_work=1
          IF(strngs(2) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(3),80)
            INQUIRE(FILE=strngs(3),EXIST=exist)
            if(strngs(3)(1:3).eq."DBG")  write_sge_work=2
            IF(exist) THEN
              CALL openf(kwork,strngs(3),'FORMATTED','OLD',0)
            ELSE
              CALL openf(kwork,strngs(3),'FORMATTED','NEW',0)
            END IF
          ELSE
            errmsg=err_args(3) //'2'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
          END IF
        ELSE if(nword.eq.4) THEN
          write_sge_work=1
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) target_replica
          IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            if(strngs(4)(1:3).eq."DBG")  write_sge_work=2
            IF(exist) THEN
              CALL openf(kwork,strngs(4),'FORMATTED','OLD',0)
            ELSE
              CALL openf(kwork,strngs(4),'FORMATTED','NEW',0)
            END IF
          ELSE
            errmsg=err_args(3) //'3'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
          END IF
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'2 or 3'
          CALL xerror(errmsg,80,1,30)
        END IF

      ELSEIF(strngs(1).EQ. 'ALWAYS_ACCEPT_PRINTF' ) THEN
         if(nword.eq.2) THEN 
            print_always=1
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nfreqsge
         else if(nword.eq.3) THEN 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nfreqsge
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) print_always
         else
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1 or 2'
          CALL xerror(errmsg,80,1,30)
        END IF

         
!---- Command BLANK LINE ----------------------------------------

    elseif ( strngs(1) == ' ' ) then
       goto 100

!---- Begininning of next ENVIRONMENT ---------------------------
    elseif ( strngs(1)(1:1) == '&' .AND. strngs(1) /= '&END' ) then
       errmsg = err_unr(1) // strngs(1)(1:8) // err_end
       call xerror(errmsg,80,1,30)
       nsevere = nsevere + 1
       goto 600

!---- Command &END -----------------------------------------------

    elseif ( strngs(1) == '&END' ) then
       goto 600
    else
       errmsg = err_unr(1) // strngs(1)(1:8) // err_end
       call xerror(errmsg,80,1,30)
       nsevere = nsevere + 1
    endif

    goto 100

600 continue

!------------------------------------------------------------
!     Environment parser ends here 
!------------------------------------------------------------

!--   Syntax errors: abort without verifying input 
    if ( nsevere > 0 .and. nsevere < 99 ) then 
       call int_str(nsevere,fmt,j)
       errmsg = fmt(1:j) //' ERRORS WHILE EXECUTING READ_SGE'
       call xerror(errmsg,80,1,2)
#ifdef _MPI_
       call mpi_barrier(mpi_comm_world,ierr) 
       call mpi_finalize(ierr)
#endif
       stop
    elseif ( nsevere > 99 ) then 
       errmsg = 'MORE THAN 99 ERRORS WHILE EXECUTING READ_SGE'
       call xerror(errmsg,80,1,2)
#ifdef _MPI_
       call mpi_barrier(mpi_comm_world,ierr) 
       call mpi_finalize(ierr)
#endif
       stop
    endif
    if ( nwarning > 0 .and. nwarning < 99 ) then 
       j = 0
       call int_str(nwarning,fmt,j)
       errmsg = fmt(1:j)//' WARNINGS WHILE EXECUTING READ_SGE'
       call xerror(errmsg,80,1,1)
    elseif ( nwarning > 99 ) then 
       errmsg = 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_SGE'
       call xerror(errmsg,80,1,1)
    endif

    return

! Errors were found

20  continue
    iret = 1
    errmsg = 'Internal reading error: wrong format? TAB character? '
    call xerror(errmsg,80,1,2)
#ifdef _MPI_
    call mpi_barrier(mpi_comm_world,ierr) 
    call mpi_finalize(ierr)
#endif
    stop

  end subroutine read_sge


!========================  ACCEPTANCE  CRITERIA  =============================================
! Application of acceptance criteria for up and down transitions between neighboring ensembles
! TSTEP = 1 means trial move for swapping from ensemble i to ensemble i+1 (e.g. increasing of the virtual temperature)
! TSTEP = -1 means trial move for swapping from ensemble i to ensemble i-1 (e.g. decreasing of the virtual temperature)
  subroutine temp_up_down ( bond_added, bend_added, tors_added, force_cost, force_ang,    &
       force_died, eqdist, eqang, eqdied, pbo, pbe, pto, da, tstep, t, radfact, lstretch, &
       lbend, litor, nbonds_added, nbends_added, nitors_added, iproc )

    use rem, only: para_index, rem_pot, rem_mat, rem_factor
    use unit, only: gascon,kwork,kprint

    implicit none

    integer, intent(in) :: tstep
    integer :: accept, i, iproc, idx
    real(8) :: delta
    real(8), external :: ranf
    real(8), intent(in) :: da, t, radfact, bond_added(*), bend_added(*), tors_added(*)
    real(8), intent(in) :: force_cost(*), force_ang(*), force_died(*)
    integer, intent(in) :: lstretch, lbend, litor, nbonds_added, nbends_added, nitors_added
    real(8) :: eqdist(*), eqang(*), eqdied(*), pbo(*), pbe(*), pto(*)

! Update the variables to calculate the acceptance ratio
    if ( tstep == 1 ) then
       tot_up(para_index) = tot_up(para_index) + 1.d0
    else
       tot_down(para_index) = tot_down(para_index) + 1.d0
    endif

    accept = 0
    delta = 0.d0

! Part related to MultiWindow Tempering on geometrical parameters (bonds and/or bends and/or torsions)
    if ( sim_tempering_sge ) then
       do i = 1, nbonds_added_st
          idx = idx_bonds(i)
          delta = delta + force_cost(idx) *                                &
               ( ( bond_added(idx) - sge_eqdist(para_index+tstep,i) )**2 - &
                 ( bond_added(idx) - sge_eqdist(para_index,      i) )**2 )
       enddo
       do i = 1, nbends_added_st
          idx = idx_bends(i)
          delta = delta + force_ang(idx) *                                        &
               ( ( bend_added(idx) - sge_eqang(para_index+tstep,i)*radfact )**2 - &
                 ( bend_added(idx) - sge_eqang(para_index,      i)*radfact )**2 )
       enddo
       do i = 1, nitors_added_st
          idx = idx_tors(i)
          delta = delta + force_died(idx) *                                        &
               ( ( tors_added(idx) - sge_eqdied(para_index+tstep,i)*radfact )**2 - &
                 ( tors_added(idx) - sge_eqdied(para_index,      i)*radfact )**2 )
       enddo
       delta = delta * 4.184
! Part related to SOLUTE TEMPERING (change of potential)
    else
       do i = 1, 3
          delta = delta + ( rem_mat(para_index+tstep,i) - rem_mat(para_index,i) ) * rem_pot(i)
       enddo
    endif

! Sum up work on a cycle between target_replicas
    if(always_accept) THEN 
       delta = delta / ( gascon * t * 0.001d0 ) 
       wtot= wtot + delta
       if(write_sge_work.eq.2) write(kwork,110) time_fs,para_index,para_index+tstep,wtot,delta
110    FORMAT(f10.1,2i5,2f12.5) 
       if(tstep == -1 .and. para_index+tstep == target_replica)  THEN  ! this is the end of the cycle
!         compute Boltz oracle. N.B if w is negative move is accepted.
          if(dexp(-wtot).lt.ranf()) THEN
             ncycle_reject=ncycle_reject+1
             use_saved=.true.
          ELSE
             go_save=.true.
          END IF
!         write(kprint,108) wtot,ncycle,ncycle_reject,use_saved,go_save
108       format(5x," SGE cycle ended ----->  Work/RT:",f10.3, " SGE-cycles: ",i12, " ( rejected ):",&
                i12," Use saved:",L1, " Go save:",L1)
          ncycle=ncycle+1
          if(write_sge_work.eq.1) write(kwork,108) wtot,ncycle,ncycle_reject,use_saved,go_save
          wtot=0.d0  ! zeroes the work 
       else
          use_saved=.false.
          go_save=.false.
       end if
    ELSE
       delta = delta / ( gascon * t * 0.001d0 ) - tstep * da
    END IF
    IF ( delta < 0. ) then
       accept = 1
    else
       if ( exp(-delta) > ranf() ) accept = 1
    endif
    if ( accept == 1 .or. always_accept) then

! Update the variables to calculate the acceptance ratio
       if ( tstep == 1 ) then
          acc_up(para_index) = acc_up(para_index) + 1.d0
       else
          acc_down(para_index) = acc_down(para_index) + 1.d0
       endif

! Updating ensemble index
       para_index = para_index + tstep

! UPDATE SCALING FACTORS
! Part related to MultiWindow Tempering on geometrical parameters
       if ( sim_tempering_sge ) then
          do i = 1, nbonds_added_st
             idx = idx_bonds(i)
             eqdist(idx) = sge_eqdist(para_index,i)
             pbo(lstretch-nbonds_added+idx) = eqdist(idx)
          enddo
          do i = 1, nbends_added_st
             idx = idx_bends(i)
             eqang(idx) = sge_eqang(para_index,i)
             pbe(lbend-nbends_added+idx) = eqang(idx) * radfact
          enddo
          do i = 1, nitors_added_st
             idx = idx_tors(i)
             eqdied(idx) = sge_eqdied(para_index,i)
             pto(litor-nitors_added+idx) = eqdied(idx) * radfact
          enddo
       else
! Part related to SOLUTE TEMPERING
          do i = 1, 3
             rem_factor(i) = rem_mat(para_index,i)
          enddo
       endif

    endif
    if(write_sge_work.eq.2) write(777,110) time_fs,para_index,para_index+tstep,wtot,delta

    return
  end subroutine temp_up_down



!==============================================================================
! Allocation of arrays useful in the Serial Generalized Ensemble algorithm
  subroutine alloc1 ( iproc )

    implicit none

    integer :: iproc
#ifdef _MPI_
    include 'mpif.h'
    integer :: ierr, status(mpi_status_size)
#endif

    allocate( sge_up(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable SGE_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    sge_up(:) = 0

    allocate( sge_up_jar(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable SGE_UP_JAR(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    sge_up_jar(:) = 0

    current_dim_up = 2
    allocate( w_up(current_dim_up,nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable W_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    w_up(:,:) = 0.

    current_dim_up_jar = 2
    allocate( w_up_jar(current_dim_up_jar,nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable W_UP_JAR(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    w_up_jar(:,:) = 0.

    allocate( sge_down(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable SGE_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    sge_down(:) = 0

    allocate( sge_down_jar(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable SGE_DOWN_JAR(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    sge_down_jar(:) = 0

    current_dim_down = 2
    allocate( w_down(current_dim_down,nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable W_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    w_down(:,:) = 0.

    current_dim_down_jar = 2
    allocate( w_down_jar(current_dim_down_jar,nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable W_DOWN_JAR(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    w_down_jar(:,:) = 0.

    allocate( estimate_df_up(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable ESTIMATE_DF_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    estimate_df_up(:) = 0

    allocate( estimate_df_down(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable ESTIMATE_DF_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    estimate_df_down(:) = 0

    allocate( df_var(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_VAR(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_var(:) = 0.d0

    allocate( df_down(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_down(:) = 0.d0

    allocate( df_up(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_up(:) = 0.d0

    allocate( nconst_down(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NCONST_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    nconst_down(:) = 0

    allocate( nconst_up(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NCONST_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    nconst_up(:) = 0

    allocate( nw_down_tot_old(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NW_DOWN_TOT_OLD(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    nw_down_tot_old(:) = 0

    allocate( nw_up_tot_old(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NW_UP_TOT_OLD(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    nw_up_tot_old(:) = 0

    allocate( acc_up(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable ACC_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    acc_up(:) = 0.

    allocate( acc_down(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable ACC_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    acc_down(:) = 0.

    allocate( tot_up(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable TOT_UP(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    tot_up(:) = 0.

    allocate( tot_down(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable TOT_DOWN(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    tot_down(:) = 0.

    allocate( histog(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable HISTOG(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    histog(:) = 0.

    allocate( histog_t(nstates), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable HISTOG_T(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    histog_t(:) = 0.

    allocate( df_fixed(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_FIXED(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_fixed(:) = 0.d0

    return
  end subroutine alloc1

!=================================================================
  subroutine alloc2 ( iproc )

    implicit none

    integer :: iproc
#ifdef _MPI_
    include 'mpif.h'
    integer :: ierr, status(mpi_status_size)
#endif

    allocate( estimate_df(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable ESTIMATE_DF(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    estimate_df(:) = 0

    allocate( dfi(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DFI(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    dfi(:) = 0.d0

    allocate( df_av(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_AV(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_av(:) = 0.d0

    allocate( ndf_av(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NDF_AV(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    ndf_av(:) = 0.d0

    allocate( ndf_av_fly(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable NDF_AV_FLY(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    ndf_av_fly(:) = 0.d0

    allocate( df(nstates-1), stat=ierror )
    if ( ierror /= 0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df(:) = 0.d0

    allocate( dfi_inst(nstates-1,sge_update), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DFI_INST(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    dfi_inst(:,:) = 0.d0

    allocate( df_var_inst(nstates-1,sge_update), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable DF_VAR_INST(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    df_var_inst(:,:) = 0.d0

    allocate( cnt_upd(nstates-1), stat=ierror )
    if ( ierror/=0 ) then
       write(*,*) 'SEVERE ERROR!!!'
       write(*,*) 'variable CNT_UPD(): allocation failed in processor', iproc
       write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
       call mpi_abort ( mpi_comm_world, ierr )
#endif
       stop
    endif
    cnt_upd(:) = 0

    return
  end subroutine alloc2

!=================================================================
  subroutine alloc3 ( iproc )

    implicit none

    integer :: iproc
#ifdef _MPI_
    include 'mpif.h'
    integer :: ierr, status(mpi_status_size)
#endif

    if ( sge_bonds ) then
       allocate(sge_eqdist(nstates,nbonds_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable SGE_EQDIST(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       sge_eqdist(:,:) = 0.d0

       allocate(idx_bonds(nbonds_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable IDX_BONDS(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       idx_bonds(:) = 0.d0
    endif

    if ( sge_bends ) then
       allocate(sge_eqang(nstates,nbends_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable SGE_EQANG(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       sge_eqang(:,:) = 0.d0

       allocate(idx_bends(nbends_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable IDX_BENDS(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       idx_bends(:) = 0.d0
    endif

    if ( sge_tors ) then
       allocate(sge_eqdied(nstates,nitors_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable SGE_EQDIED(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       sge_eqdied(:,:) = 0.d0

       allocate(idx_tors(nitors_added_st), stat=ierror)
       if ( ierror /= 0 ) then
          write(*,*) 'SEVERE ERROR!!!'
          write(*,*) 'variable IDX_TORS(): allocation failed in processor', iproc
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
       idx_tors(:) = 0.d0
    endif

    return
  end subroutine alloc3


!=============================  MAKES INITIALIZATIONS  ========================================
  subroutine sge_init( flag, iproc, nproc, nbonds_added, nbends_added, nitors_added, ntap, ltors, int14p, ltor, int14, &
       xp0, yp0, zp0, radfact, force_cost, eqdist, eqdist1, force_ang, eqang, eqang1, force_died, eqdied, eqdied1,     &
       atom_b1, atom_b2, atom_be1, atom_be2, atom_be3,atom_it1, atom_it2, atom_it3, atom_it4, strbonds, strbends,      &
       strtors, remltor, remint14, rematom )

    use rem, only: rem_group, rem_groups, rem_mat, rem_factor, rem_segment, rem_factor_max, para_index, cdist

    implicit none

    real(8), intent(in) :: xp0(*), yp0(*), zp0(*)
    real(8), intent(in) :: radfact
    integer :: nbonds_added_tmp, nbends_added_tmp, nitors_added_tmp
    integer :: nbonds_added, nbends_added, nitors_added, nc
    integer :: atom1, atom2, atom3, atom4
    integer :: atom_b1(*), atom_b2(*)
    integer :: atom_be1(*), atom_be2(*), atom_be3(*)
    integer :: atom_it1(*), atom_it2(*), atom_it3(*), atom_it4(*)
    real(8) :: force_cost(*), eqdist(*), eqdist1(*)
    real(8) :: force_ang(*), eqang(*), eqang1(*)
    real(8) :: force_died(*), eqdied(*), eqdied1(*)
    real(8) :: fc, emin
    integer :: strbonds(*), strbends(*), strtors(*)
    integer :: flag, iproc, nproc, ntap, ltors, int14p
    integer :: ltor(4,*), int14(2,*)
    logical :: remltor(*), remint14(*), rematom(*)
    integer :: i, i1, j, itmp, idx, jmin
    integer :: nerrors, ierr
    real(8) :: dum1, dum2, dum3, delta_coord, xxx
    real(8) :: buffer1(0:100), buffer2(0:100), buffer3(0:100)
    character(27) :: err_unr(4)
    character(15) :: err_fnf
    character(80) :: errmsg
    logical :: exist
    real(8), external :: ranf

! Random number generator is initialized
    do i = 0, iproc
       xxx = ranf()
    enddo

    err_unr(1)='Unrecognized command  ---> '
    err_unr(2)='Unrecognized subcommand -> '
    err_unr(3)='Unrecognized keyword ----> '
    err_unr(4)='UNSUPPORTED  COMMAND ----> '
    err_fnf = ' file not found' 
    nerrors = 0

! Allocation of arrays useful in Serial Generalized Ensemble algorithm
    call alloc1 ( iproc )
    if ( flag == 0 ) call alloc2 ( iproc )

! Set nt_points_w: used in sge_df_update subroutine
    nt_points_w = int(sge_tsave_pmf / sge_tsave_work)

! sge_df_opt = 0: use fixed free energies
! sge_df_opt = 1: update free energies (using BAR-SGE) during the run -- (JCTC 2010, 6, 1935-1950)
    if ( sge_df_opt == 0 ) then
       open( unit=999, file=df_file, form='formatted' )
       do i = 1, nstates - 1
          read(999,*,end=1005) df_fixed(i)
       enddo
       close(999)
       goto 1006
1005   errmsg = 'End of file encountered in ' // df_file
       call xerror(errmsg,80,1,30)
       errmsg = 'The number of free energy differences in the file'
       call xerror(errmsg,80,1,25)
       errmsg = df_file
       call xerror(errmsg,80,1,25)
       errmsg = 'must be equal to the number of states/ensembles minus 1'
       call xerror(errmsg,80,1,25)
       nerrors = nerrors + 1
       goto 1600
1006   continue
    endif

! cdist = 0: restart from an old run
! cdist = 1: new run with a linear progression for reaction coordinate
! cdist = 2: new run with a custom progression

    if ( cdist == 0 .or. cdist == 2 ) then
#ifdef _MPI_
       if ( cdist == 0 ) then
          inquire( file='SGE.set', exist=exist )
       else
          inquire( file='../SGE.set', exist=exist )
       endif
#else
       inquire( file='SGE.set', exist=exist )
#endif
       if ( .not.exist ) then
#ifdef _MPI_
          if ( cdist == 0 ) then
             errmsg = 'Cannot find SGE.set for restart'
          else
             errmsg = 'Cannot find ../SGE.set for restart'
          endif
#else
          errmsg = 'Cannot find SGE.set for restart'
#endif
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
          goto 1600
       endif
#ifdef _MPI_
       if ( cdist == 0 ) then
          open( unit=999, file='SGE.set', form='formatted' )
       else
          open( unit=999, file='../SGE.set', form='formatted' )
       endif
#else
       open( unit=999, file='SGE.set', form='formatted' )
#endif
    endif


    if ( sim_tempering_sge ) then

!-----------------------------------------------------------------------------
! Set envirnment for MultiWindow Tempering simulation in the space of collective coordinates
       sge_bonds = .false.
       sge_bends = .false.
       sge_tors = .false.

       nbonds_added_st = 0
       do i = 1, nbonds_added
          if ( strbonds(i).EQ.1 ) then
             sge_bonds = .true.
             nbonds_added_st = nbonds_added_st + 1
          endif
       enddo

       nbends_added_st = 0
       do i = 1, nbends_added
          if ( strbends(i).EQ.1 ) then
             sge_bends = .true.
             nbends_added_st = nbends_added_st + 1
          endif
       enddo

       nitors_added_st = 0
       do i = 1, nitors_added
          if ( strtors(i) .EQ. 1) then
             sge_tors = .true.
             nitors_added_st = nitors_added_st + 1
          endif
       enddo

       if ( .not.sge_bonds .and. .not.sge_bends .and. .not.sge_tors ) then
          errmsg = 'In MultiWindow Tempering simulations,'
          call xerror(errmsg,80,1,30)
          errmsg = 'at least one bonding/bending/torsion must be'
          call xerror(errmsg,80,1,26)
          errmsg = 'defined using ADD_STR_X (X=BONDS,BENDS,TORS)'
          call xerror(errmsg,80,1,26)
          nerrors = nerrors + 1
          goto 1600
       endif

       call alloc3 ( iproc )

! restart an old run or restart from custom file
       if ( cdist == 0 .or. cdist == 2 ) then

          read(999,*)
          read(999,*) nbonds_added_tmp,nbends_added_tmp,nitors_added_tmp

          if ( nbonds_added_tmp /= nbonds_added_st ) then
             if ( cdist == 0 ) then
                errmsg = 'In SGE.set the number of added bondings'
                call xerror(errmsg,80,1,30)
                errmsg = '(see first field in first line of SGE.set)'
             else
#ifdef _MPI_
                errmsg = 'In ../SGE.set the number of added bondings'
#else
                errmsg = 'In SGE.set the number of added bondings'
#endif
                call xerror(errmsg,80,1,30)
#ifdef _MPI_
                errmsg = '(see first field in first line of ../SGE.set)'
#else
                errmsg = '(see first field in first line of SGE.set)'
#endif
             endif
             call xerror(errmsg,80,1,26)
             errmsg = 'does not correspond to that defined in'
             call xerror(errmsg,80,1,26)
             errmsg='ADD_STR_BONDS(&POTENTIAL)'
             call xerror(errmsg,80,1,26)
             nerrors = nerrors + 1
             goto 1600
          endif

          if ( nbends_added_tmp /= nbends_added_st ) then
             if ( cdist == 0 ) then
                errmsg = 'In SGE.set the number of added bendings'
                call xerror(errmsg,80,1,30)
                errmsg = '(see second field in first line of SGE.set)'
             else
#ifdef _MPI_
                errmsg = 'In ../SGE.set the number of added bendings'
#else
                errmsg = 'In SGE.set the number of added bendings'
#endif
                call xerror(errmsg,80,1,30)
#ifdef _MPI_
                errmsg = '(see second field in first line of ../SGE.set)'
#else
                errmsg = '(see second field in first line of SGE.set)'
#endif
             endif
             call xerror(errmsg,80,1,26)
             errmsg = ' does not correspond to that defined in'
             call xerror(errmsg,80,1,26)
             errmsg='ADD_STR_BENDS(&POTENTIAL)'
             call xerror(errmsg,80,1,26)
             nerrors = nerrors + 1
             goto 1600
          endif

          if ( nitors_added_tmp /= nitors_added_st ) then
             if ( cdist == 0 ) then
                errmsg = 'In SGE.set the number of added torsions'
                call xerror(errmsg,80,1,30)
                errmsg = '(see third field in first line of SGE.set)'
             else
#ifdef _MPI_
                errmsg = 'In ../SGE.set the number of added torsions'
#else
                errmsg = 'In SGE.set the number of added torsions'
#endif
                call xerror(errmsg,80,1,30)
#ifdef _MPI_
                errmsg = '(see third field in first line of ../SGE.set)'
#else
                errmsg = '(see third field in first line of SGE.set)'
#endif
             endif
             call xerror(errmsg,80,1,26)
             errmsg = 'does not correspond to that defined in'
             call xerror(errmsg,80,1,26)
             errmsg = 'ADD_STR_TORS(&POTENTIAL)'
             call xerror(errmsg,80,1,26)
             nerrors = nerrors + 1
             goto 1600
          endif

          do j = 1, nbonds_added_st
             read(999,*)
             read(999,*) atom1, atom2, fc, idx
             if ( strbonds(idx) .EQ. 1) then
                if ( atom1 /= atom_b1(idx) .or. atom2 /= atom_b2(idx) .or. dabs(fc-force_cost(idx)) > 1.d-8 ) then
#ifdef _MPI_
                   if ( cdist == 0 ) then
                      errmsg = 'In SGE.set the bonding atoms or forces'
                   else
                      errmsg = 'In ../SGE.set the bonding atoms or forces'
                   endif
#else
                   errmsg = 'In SGE.set the bonding atoms or forces'
#endif
                   call xerror(errmsg,80,1,30)
                   errmsg = 'do not match the input from'
                   call xerror(errmsg,80,1,26)
                   errmsg = 'ADD_STR_BONDS(&POTENTIAL)'
                   call xerror(errmsg,80,1,26)
                   nerrors = nerrors + 1
                   goto 1600
                endif
             else
#ifdef _MPI_
                if ( cdist == 0 ) then
                   errmsg = 'In SGE.set the bonding-index does not'
                else
                   errmsg = 'In ../SGE.set the bonding-index does not'
                endif
#else
                errmsg = 'In SGE.set the bonding-index does not'
#endif
                call xerror(errmsg,80,1,30)
                errmsg = 'correspond to that resulting from'
                call xerror(errmsg,80,1,26)
                errmsg = 'ADD_STR_BONDS(&POTENTIAL)'
                call xerror(errmsg,80,1,26)
                nerrors = nerrors + 1
                goto 1600
             endif
             idx_bonds(j) = idx
             read(999,*)
             do i = 1, nstates
                read(999,*,end=2003) itmp, sge_eqdist(i,j)
             enddo
          enddo

          do j = 1, nbends_added_st
             read(999,*)
             read(999,*) atom1, atom2, atom3, fc, idx
             if ( strbends(idx).EQ. 1 ) then
                if ( atom1 /= atom_be1(idx) .or. atom2 /= atom_be2(idx) .or. &
                     atom3 /= atom_be3(idx) .or. dabs(fc-force_ang(idx)) > 1.d-8 ) then
#ifdef _MPI_
                   if ( cdist == 0 ) then
                      errmsg = 'In SGE.set the bending atoms or forces'
                   else
                      errmsg = 'In ../SGE.set the bending atoms or forces'
                   endif
#else
                   errmsg = 'In SGE.set the bending atoms or forces'
#endif
                   call xerror(errmsg,80,1,30)
                   errmsg = 'do not match the input from'
                   call xerror(errmsg,80,1,26)
                   errmsg = 'ADD_STR_BENDS(&POTENTIAL)'
                   call xerror(errmsg,80,1,26)
                   nerrors = nerrors + 1
                   goto 1600
                endif
             else
#ifdef _MPI_
                if ( cdist == 0 ) then
                   errmsg = 'In SGE.set the bending-index does not'
                else
                   errmsg = 'In ../SGE.set the bending-index does not'
                endif
#else
                errmsg = 'In SGE.set the bending-index does not'
#endif
                call xerror(errmsg,80,1,30)
                errmsg = 'correspond to that resulting from'
                call xerror(errmsg,80,1,26)
                errmsg = 'ADD_STR_BENDS(&POTENTIAL)'
                call xerror(errmsg,80,1,26)
                nerrors = nerrors + 1
                goto 1600
             endif
             idx_bends(j) = idx
             read(999,*)
             do i = 1, nstates
                read(999,*,end=2003) itmp, sge_eqang(i,j)
             enddo
          enddo

          do j = 1, nitors_added_st
             read(999,*)
             read(999,*) atom1, atom2, atom3, atom4, fc, idx
             if ( strtors(idx) .EQ. 1) then
                if ( atom1 /= atom_it1(idx) .or. atom2 /= atom_it2(idx) .or. &
                     atom3 /= atom_it3(idx) .or. atom4 /= atom_it4(idx) .or. &
                     dabs(fc-force_died(idx)) > 1.d-8 ) then
#ifdef _MPI_
                   if ( cdist == 0 ) then
                      errmsg = 'In SGE.set the torsion atoms or forces'
                   else
                      errmsg = 'In ../SGE.set the torsion atoms or forces'
                   endif
#else
                   errmsg = 'In SGE.set the torsion atoms or forces'
#endif
                   call xerror(errmsg,80,1,30)
                   errmsg = 'do not match the input from'
                   call xerror(errmsg,80,1,26)
                   errmsg = 'ADD_STR_TORS(&POTENTIAL)'
                   call xerror(errmsg,80,1,26)
                   nerrors = nerrors + 1
                   goto 1600
                endif
             else
#ifdef _MPI_
                if ( cdist == 0 ) then
                   errmsg = 'In SGE.set the torsion-index does not'
                else
                   errmsg = 'In ../SGE.set the torsion-index does not'
                endif
#else
                errmsg = 'In SGE.set the torsion-index does not'
#endif
                call xerror(errmsg,80,1,30)
                errmsg = 'correspond to that resulting from'
                call xerror(errmsg,80,1,26)
                errmsg = 'ADD_STR_TORS(&POTENTIAL)'
                call xerror(errmsg,80,1,26)
                nerrors = nerrors + 1
                goto 1600
             endif
             idx_tors(j) = idx
             read(999,*)
             do i = 1, nstates
                read(999,*,end=2003) itmp, sge_eqdied(i,j)
             enddo
          enddo

          close(999)

          goto 2004

#ifdef _MPI_
2003      if ( cdist == 0 ) then
             errmsg = 'End of file encountered in SGE.set'
          else
             errmsg = 'End of file encountered in ../SGE.set'
          endif
#else
2003      errmsg = 'End of file encountered in SGE.set'
#endif
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
          goto 1600

2004      continue

       elseif ( cdist == 1 ) then
! use a linear progression for the reaction coordinate

          if ( sge_bonds ) then
             nc = 0
             do i = 1, nbonds_added
                if ( strbonds(i).EQ. 1 ) then
                   nc = nc + 1
                   idx_bonds(nc) = i
                   delta_coord = (eqdist1(i) - eqdist(i)) / (nstates-1)
                   sge_eqdist(1,nc) = eqdist(i)
                   do j = 2, nstates
                      sge_eqdist(j,nc) = sge_eqdist(j-1,nc) + delta_coord
                   enddo
                endif
             enddo
          endif

          if ( sge_bends ) then
             nc = 0
             do i = 1, nbends_added
                if ( strbends(i) .EQ. 1) then
                   nc = nc + 1
                   idx_bends(nc) = i
                   delta_coord = (eqang1(i) - eqang(i)) / (nstates-1)
                   sge_eqang(1,nc) = eqang(i)
                   do j = 2, nstates
                      sge_eqang(j,nc) = sge_eqang(j-1,nc) + delta_coord
                   enddo
                endif
             enddo
          endif

          if ( sge_tors ) then
             nc = 0
             do i = 1, nitors_added
                if ( strtors(i).EQ. 1 ) then
                   nc = nc + 1
                   idx_tors(nc) = i
                   delta_coord = (eqdied1(i) - eqdied(i)) / (nstates-1)
                   sge_eqdied(1,nc) = eqdied(i)
                   do j = 2, nstates
                      sge_eqdied(j,nc) = sge_eqdied(j-1,nc) + delta_coord
                   enddo
                endif
             enddo
          endif

       else

          nerrors = nerrors + 1
          errmsg = err_unr(3) // 'in SETUP(&SGE): must be 0, 1 or 2'
          call xerror(errmsg,80,1,30)
          goto 1600

       endif

! Multi-Window Tempering: assignment of the replicas to the ensemble with lowest initial energy (only for cold runs)
       if ( cdist /= 0 ) then
          call init_bond_bend_tors ( jmin, emin, radfact, xp0, yp0, zp0, force_cost, force_ang, force_died, eqdied, &
               atom_b1, atom_b2, atom_be1, atom_be2, atom_be3, atom_it1, atom_it2, atom_it3, atom_it4 )
          para_index = jmin
       endif

       do i = 1, nbonds_added_st
          eqdist(idx_bonds(i)) = sge_eqdist(para_index,i)
       enddo

       do i = 1, nbends_added_st
          eqang(idx_bends(i)) = sge_eqang(para_index,i)
       enddo

       do i = 1, nitors_added_st
          eqdied(idx_tors(i)) = sge_eqdied(para_index,i)
       enddo

#ifdef _MPI_
       if ( cdist > 0 ) then
#else
       if ( cdist == 1 ) then
#endif
          open( unit=999, file='SGE.set', form='formatted' )
          write(999,'(3(a11,2x))') 'Added-bonds', 'Added-bends', ' Added-tors'
          write(999,'(3(i11,2x))') nbonds_added_st, nbends_added_st, nitors_added_st

          if ( nbonds_added_st > 0 ) then
             do j = 1, nbonds_added_st
                idx = idx_bonds(j)
                write(999,'(4(a14,4x))') '        Atom-1', '        Atom-2', 'Force-constant', '         Index'
                write(999,1500) atom_b1(idx), atom_b2(idx), force_cost(idx), idx
1500            format(2(i14,4x),f14.5,4x,i14)
                write(999,'(a10,10x,a22)') 'n.Ensemble', 'Coordinate (Eq. value)'
                do i = 1, nstates
                   write(999,'(i10,2x,g30.15)') i, sge_eqdist(i,j)
                enddo
             enddo
          endif

          if ( nbends_added_st > 0 ) then
             do j = 1, nbends_added_st
                idx = idx_bends(j)
                write(999,'(5(a14,4x))') '        Atom-1', '        Atom-2', '        Atom-3', 'Force-constant', '         Index'
                write(999,1501) atom_be1(idx), atom_be2(idx), atom_be3(idx), force_ang(idx), idx
1501            format(3(i14,4x),f14.5,4x,i14)
                write(999,'(a10,10x,a22)') 'n.Ensemble', 'Coordinate (Eq. value)'
                do i = 1, nstates
                   write(999,'(i10,2x,g30.15)') i, sge_eqang(i,j)
                enddo
             enddo
          endif

          if ( nitors_added_st > 0 ) then
             do j = 1, nitors_added_st
                idx = idx_tors(j)
                write(999,'(6(a14,4x))') &
                     '        Atom-1', '        Atom-2', '        Atom-3', '        Atom-4', 'Force-constant', '         Index'
                write(999,1502) atom_it1(idx), atom_it2(idx), atom_it3(idx), atom_it4(idx), force_died(idx), idx
1502            format(4(i14,4x),f14.5,4x,i14)
                write(999,'(a10,10x,a22)') 'n.Ensemble', 'Coordinate (Eq. value)'
                do i = 1, nstates
                   write(999,'(i10,2x,g30.15)') i, sge_eqdied(i,j)
                enddo
             enddo
          endif

          close(999)
       endif

    else

!---------------------------------------------------------------------------------
! SET ENVIRNMENT FOR SGE HAMILTONIAN SIMULATION (solute tempering)

! cdist = 0: restart from an old run
! cdist = 1: new run with a geometric progression for temperatures
! cdist = 2: new run with a custom progression from ../SGE.set file

       if ( cdist == 0 .or. cdist == 2 ) then

          read(999,*)
          read(999,*)
          do i = 1, nstates
             read(999,*,end=1001) itmp, ( rem_mat(i,j), j=1, 3 )
          enddo
          close(999)

          goto 1002
1001      errmsg = 'End of file encountered in SGE.set'
          call xerror(errmsg,80,1,30)
          nerrors = nerrors + 1
          goto 1600

1002      continue

       elseif ( cdist == 1 ) then ! use a geometric progression

          do j = 1, nstates
             do i = 1, 3
                rem_mat(j,i) = rem_factor_max(i)**( dfloat(j-1) / (nstates - 1.d0) )
             enddo
          enddo

       else
          
          nerrors = nerrors + 1
          errmsg = err_unr(3) // 'in SETUP(&SGE): must be 0, 1 or 2'
          call xerror(errmsg,80,1,30)
          goto 1600

       endif

! assign to each process a replica index and a scaling factor
       if ( cdist == 0 ) then

          do i = 1, 3
             rem_factor(i) = rem_mat(para_index,i)
          enddo
!          write(*,8099) para_index, ( rem_factor(i), i=1, 3 )
!8099      format( '-- FROM RESTART FILE:', /, '-- REPLICA INDEX',I10,/,'-- SCALING_FACTOR',3F10.4/) 

       else

          para_index = 0
          do i = 0, nproc-1
             para_index = para_index + 1
             if ( para_index > nstates ) para_index = 1
             if ( i == iproc ) goto 150
          enddo

150       do i = 1, 3
             rem_factor(i) = rem_mat(para_index,i)
          enddo
#ifdef _MPI_
          if ( cdist > 0 ) then
#else
          if ( cdist == 1 ) then
#endif
             open( unit=999, file='SGE.set', form='formatted' )
             write(999,'(a)') 'Scaling factors:'
             write(999,'(a)') 'n.Ensemble      Bending+Bonding               Torsions+1-4'//'                  Non-Bonded'
             do i = 1, nstates
                write(999,'(i6,1x,3g30.15)') i, ( rem_mat(i,j), j=1, 3 )
             enddo
             close(999)
          endif

       endif

! **************************  prepare a solute tempering run  ****************************

       if ( rem_segment ) then
! build three logical arrays, remltor and remint14 and rematom
          do i = 1,int14p
             remint14(i) = .false.
          enddo
          do i = 1,ltors
             remltor(i) = .false.
             do j = 1,rem_groups
                do i1 = 1,4
! if the atoms of the torsion overlap with the subsets of atoms, then 
                   if ((ltor(i1,i) >= rem_group(1,j)) .and. (ltor(i1,i) <= rem_group(2,j))) remltor(i) = .true.
                enddo
             enddo
             if (remltor(i)) then 
! the 1-4 interaction will be scaled too
! find the corresponding 1-4 interaction
                do j = 1,int14p
                   if ((int14(1,j) == ltor(1,i)) .and. (int14(2 ,j) == ltor(4,i))) then
                      remint14(j) = .true.
                   endif
                enddo
             endif
          enddo
          do i = 1,ntap
             rematom(i) = .false.
             do j = 1,rem_groups
                if ((i.ge.rem_group(1,j)).and.(i.le.rem_group(2,j))) rematom(i) = .true.
             enddo
          enddo
       endif

    endif

! if errors were found, stop
1600 if ( nerrors > 0 ) then
       errmsg= ' Errors while verifying SGE input'
       call xerror(errmsg,80,1,2)
       stop
    endif

    return 
  end subroutine sge_init


!================  CALCULATE INITIAL STRETCHING, BENDINGS and TORSIONS for MULTI WINDOW TEMPERING  ============
  subroutine init_bond_bend_tors ( jmin, emin, radfact, xp0, yp0, zp0, force_cost, force_ang, force_died, eqdied, &
       atom_b1, atom_b2, atom_be1, atom_be2, atom_be3, atom_it1, atom_it2, atom_it3, atom_it4 )

    implicit none

    integer :: i, j, idx
    real(8) :: bond_val, bend_val, tors_val, et, emin
    real(8) :: x12, y12, z12, x23, y23, z23, x21, y21, z21, x32, y32, z32, x43, y43 &
         ,z43, rs12, rs32, cb, rsq21, rsq32, rsq43, rsp21, rsp32, rsp43, cb1, cb2   &
         ,cb3, sb1, sb2, sb3, aux, coa, xr1, xr2, xr3, xr4, yr1, yr2, yr3, yr4, zr1 &
         ,zr2, zr3, zr4, xn, yn, zn
    real(8) :: eqdied(*)
    real(8), intent(in) :: radfact
    real(8), intent(in) :: xp0(*), yp0(*), zp0(*)
    real(8), intent(in) :: force_cost(*), force_ang(*), force_died(*)
    integer, intent(in) :: atom_b1(*), atom_b2(*)
    integer, intent(in) :: atom_be1(*), atom_be2(*), atom_be3(*)
    integer, intent(in) :: atom_it1(*), atom_it2(*), atom_it3(*), atom_it4(*)
    integer, intent(out) :: jmin

    do j = 1, nstates
       et = 0.d0
!--   compute bonds
       do i = 1, nbonds_added_st
          idx = idx_bonds(i)
          x12 = xp0(atom_b1(idx)) - xp0(atom_b2(idx))
          y12 = yp0(atom_b1(idx)) - yp0(atom_b2(idx))
          z12 = zp0(atom_b1(idx)) - zp0(atom_b2(idx))
          bond_val = dsqrt( x12**2 + y12**2 + z12**2 )
          et = et + force_cost(idx) * ( bond_val - sge_eqdist(j,i) )**2
!          write(*,'(a,2i4,3f10.4)') 'bond ', j, i, bond_val, sge_eqdist(j,i), et
       enddo
!--   compute bends
       do i = 1, nbends_added_st
          idx = idx_bends(i)
          x12 = xp0(atom_be1(idx)) - xp0(atom_be2(idx))
          y12 = yp0(atom_be1(idx)) - yp0(atom_be2(idx))
          z12 = zp0(atom_be1(idx)) - zp0(atom_be2(idx))
          x23 = xp0(atom_be3(idx)) - xp0(atom_be2(idx))
          y23 = yp0(atom_be3(idx)) - yp0(atom_be2(idx))
          z23 = zp0(atom_be3(idx)) - zp0(atom_be2(idx))
          rs12 = x12**2 + y12**2 + z12**2
          rs32 = x23**2 + y23**2 + z23**2
          cb = ( x12 * x23 + y12 * y23 + z12 * z23 ) / dsqrt( rs12 * rs32 )
          bend_val = DACOS(cb)
          et = et + force_ang(idx) * ( bend_val - sge_eqang(j,i) * radfact )**2
!          write(*,'(a,2i4,3f10.4)') 'bend ', j, i, bend_val / radfact, sge_eqang(j,i), et
       enddo
!--   compute torsions
       do i = 1, nitors_added_st
          idx = idx_tors(i)
          eqdied(idx) = sge_eqdied(j,i) * radfact
          xr1 = xp0(atom_it1(idx)) 
          xr2 = xp0(atom_it2(idx)) 
          xr3 = xp0(atom_it3(idx)) 
          xr4 = xp0(atom_it4(idx)) 
          yr1 = yp0(atom_it1(idx)) 
          yr2 = yp0(atom_it2(idx)) 
          yr3 = yp0(atom_it3(idx)) 
          yr4 = yp0(atom_it4(idx)) 
          zr1 = zp0(atom_it1(idx)) 
          zr2 = zp0(atom_it2(idx)) 
          zr3 = zp0(atom_it3(idx)) 
          zr4 = zp0(atom_it4(idx)) 
          x21 = xr2 - xr1
          y21 = yr2 - yr1
          z21 = zr2 - zr1
          x32 = xr3 - xr2
          y32 = yr3 - yr2
          z32 = zr3 - zr2
          x43 = xr4 - xr3
          y43 = yr4 - yr3
          z43 = zr4 - zr3
          xn = y21 * z32 - y32 * z21
          yn = z21 * x32 - z32 * x21
          zn = x21 * y32 - x32 * y21
          rsq21 = x21**2 + y21**2 + z21**2
          rsq32 = x32**2 + y32**2 + z32**2
          rsq43 = x43**2 + y43**2 + z43**2
          rsp21 = dsqrt(rsq21)
          rsp32 = dsqrt(rsq32)
          rsp43 = dsqrt(rsq43)
          cb1 = ( x21*x32 + y21*y32 + z21*z32 ) / ( rsp21*rsp32 )
          cb2 = ( x43*x32 + y43*y32 + z43*z32 ) / ( rsp43*rsp32 )
          cb3 = ( x21*x43 + y21*y43 + z21*z43 ) / ( rsp21*rsp43 )
          sb1 = dsqrt( 1.0d0 - cb1**2 )
          sb2 = dsqrt( 1.0d0 - cb2**2 )
          sb3 = dsqrt( 1.0d0 - cb3**2 )
          aux = sb1 * sb2
          coa = ( cb1 * cb2 - cb3 ) / aux
          tors_val = dacos(coa)
          aux = sign( tors_val, xn * x43 + yn * y43 + zn * z43 )
          if( aux < 0.d0 .and. dabs(aux) > 1.57079 ) then
             tors_val = 2 * dacos(-1.d0) + aux
             if( abs( eqdied(idx) - tors_val ) > 3.14159D0 ) tors_val = -2 * 3.14159D0 + tors_val
          elseif( aux < 0.d0 .and. dabs(aux) <= 1.57079 ) then
             tors_val = aux
             if( abs( eqdied(idx) - tors_val ) > 3.14159D0 ) tors_val = 2 * 3.14159D0 + tors_val
          else
             tors_val = aux
             if( abs( eqdied(idx) - tors_val ) > 3.14159D0 ) tors_val = 2 * 3.14159D0 + tors_val
             if( abs( eqdied(idx) - tors_val ) > 3.14159D0 ) tors_val = -4 * 3.14159D0 + tors_val
          endif
          et = et + force_died(idx) * ( tors_val - eqdied(idx) )**2
!          write(*,'(a,2i4,3f10.4)') 'tors ', j, i, tors_val / radfact, eqdied(idx)/ radfact, et
       enddo
       if ( j == 1 ) then
          emin = et
          jmin = j
       else
          if ( et < emin ) then
             emin = et
             jmin = j
          endif
       endif
    enddo

    return
  end subroutine init_bond_bend_tors


!=============================  PRINT PARAMETERS AT THE BEGINNING OF THE RUN  ========================
  subroutine sge_print_titles ( iproc, nproc, atom_b1, atom_b2, atom_be1, atom_be2, atom_be3, atom_it1, atom_it2, &
       atom_it3, atom_it4 )

    use cpropar, only: t, force_cost, force_ang, force_died
    use unit, only: kprint
    use rem, only: cdist, rem_factor, rem_factor_max, para_index, rem_groups, rem_segkind, rem_group

    implicit none

    integer, intent(in) :: iproc, nproc
    integer, intent(in) :: atom_b1(*), atom_b2(*)
    integer, intent(in) :: atom_be1(*), atom_be2(*), atom_be3(*)
    integer, intent(in) :: atom_it1(*), atom_it2(*), atom_it3(*), atom_it4(*)

    integer :: i, idx, tlscaled
    character(90) :: line_res1, line_res2, sge_line_res

    if ( cdist == 0 ) then
       line_res1 = 'Simulation restarted from a previous run'
       line_res2 = '     The ensemble-parameters are taken from SGE.set'
    elseif ( cdist == 1 ) then
       if ( sim_tempering_sge ) then
          line_res1 = 'A linear progression is employed for the reaction coordinate (default)'
       else
          line_res1 = 'A geometric progression is employed for the ensemble-temperatures (default)'
       endif
       line_res2 = '     The ensemble-parameters are written to the SGE.set file'
    else
       if ( sim_tempering_sge ) then
          line_res1 = 'A customized progression is employed for the reaction coordinate'
       else
          line_res1 = 'A customized progression is employed for the ensemble-temperatures'
       endif
#ifdef _MPI_
       line_res2 = '     The ensemble-parameters are taken from ../SGE.set'
#else
       line_res2 = '     The ensemble-parameters are taken from SGE.set'
#endif
    endif

    if ( sge_df_opt == 0 ) then
       sge_line_res = 'Use fixed free energies taken from file ' // df_file
    else
       sge_line_res = 'Use BAR-SGE method to update free energies'
    endif

    write(kprint,8100) nproc, nstates, iproc, para_index, sge_transition_scheme, sge_ts, sge_tsave_work, &
         sge_tsave_pmf, sge_update, sge_f, sge_print

    if ( .not. sim_tempering_sge ) then
       write(kprint,'(a,f10.2)') ' ------ Lowest temperature:              ', t
       write(kprint,'(a,3f10.2)') ' ------ Highest temperatures:            ', &
            t / rem_factor_max(1), t / rem_factor_max(2), t / rem_factor_max(3)
       write(kprint,'(a,3f10.2)') ' ------ Current temperatures:            ', &
            t / rem_factor(1), t / rem_factor(2), t / rem_factor(3)
    endif

    write(kprint,'(a,a90)') ' ------ ', line_res1
    write(kprint,'(a,a90)') ' ------ ', line_res2
    write(kprint,'(a,a90)') ' ------ ', sge_line_res

    if ( sge_df_opt == 0 ) then
       write(kprint,*) '       Fixed dimensionless free energy differences:'
       do i = 1, nstates-1
          write(kprint,'(5x,i4,1x,a3,i4,5x,a4,1x,f12.5)') i, '-->', i+1, 'Df =', df_fixed(i)
       enddo
    endif

    if ( rem_groups >= 1 ) then
       write(kprint,8001)
       tlscaled = 0
       do i = 1, rem_groups
          write(kprint,8002) i, rem_group(1,i), rem_group(2,i)
          tlscaled = tlscaled + rem_group(2,i) - rem_group(1,i) + 1
       enddo
       write(kprint,8003) tlscaled
       if ( rem_segkind == 0 ) write(kprint,8004)
       if ( rem_segkind == 1 ) write(kprint,8005)
       if ( rem_segkind == 2 ) write(kprint,8006)
    endif

#ifdef _MPI_
    if ( iproc == 0 ) then
       write(*,8101) nproc, nstates, sge_transition_scheme, sge_ts, sge_tsave_work, sge_tsave_pmf, sge_update, sge_f, sge_print

       if ( .not. sim_tempering_sge ) then
          write(*,'(a,f10.2)') ' ------ Lowest temperature:              ', t
          write(*,'(a,3f10.2)') ' ------ Highest temperatures:            ', &
               t/rem_factor_max(1), t/rem_factor_max(2), t/rem_factor_max(3)
          write(*,'(a,3f10.2)') ' ------ Current temperatures:            ', &
               t/rem_factor(1), t /rem_factor(2), t/rem_factor(3)
       endif

       write(*,'(a,a90)') ' ------ ', line_res1
       write(*,'(a,a90)') ' ------ ', line_res2
       write(*,'(a,a90)') ' ------ ', sge_line_res

       if ( sge_df_opt == 0 ) then
          write(*,*) '       Fixed dimensionless free energy differences:'
          do i = 1, nstates-1
             write(*,'(5x,i4,1x,a3,i4,5x,a4,1x,f12.5)') i, '-->', i+1, 'Df =', df_fixed(i)
          enddo
       endif

       if ( rem_groups >= 1 ) then
          write(*,8001)
          tlscaled = 0
          do i = 1, rem_groups
             write(*,8002) i, rem_group(1,i), rem_group(2,i)
             tlscaled = tlscaled + rem_group(2,i) - rem_group(1,i) + 1
          enddo
          write(*,8003) tlscaled
          if ( rem_segkind == 0 ) write(*,8004)
          if ( rem_segkind == 1 ) write(*,8005)
          if ( rem_segkind == 2 ) write(*,8006)
       endif
    endif
#endif

    if ( sim_tempering_sge ) then

       write(kprint,*)
#ifdef _MPI_
       if ( iproc == 0 ) write(*,*)
#endif

       if ( sge_bonds ) then
          write(kprint,'(a)') 'Bonding reaction-coordinates in Multi-Window Tempering'
          write(kprint,1580)
          do i = 1, nbonds_added_st
             idx = idx_bonds(i)
             write(kprint,1581) atom_b1(idx), atom_b2(idx), force_cost(idx), sge_eqdist(1,i), sge_eqdist(nstates,i)
          enddo
#ifdef _MPI_
          if ( iproc == 0 ) then
             write(*,'(a)') 'Bonding reaction-coordinates in Multi-Window Tempering'
             write(*,1580)
             do i = 1, nbonds_added_st
                idx = idx_bonds(i)
                write(*,1581) atom_b1(idx), atom_b2(idx), force_cost(idx), sge_eqdist(1,i), sge_eqdist(nstates,i)
             enddo
          endif
#endif
       endif

1580   format(5x, 'Atom1', 5x, 'Atom2', 5x, 'K (kcal/mol/Ang**2)', 1x, '           r0 (Ang)', 1x, '           r1 (Ang)' )
1581   format( 2i10, 4x, 3(8x,f12.3) )

       if ( sge_bends ) then
          write(kprint,'(a)') 'Bending reaction-coordinates in Multi-Window Tempering'
          write(kprint,1582)
          do i = 1, nbends_added_st
             idx = idx_bends(i)
             write(kprint,1583) atom_be1(idx), atom_be2(idx), atom_be3(idx), force_ang(idx), sge_eqang(1,i), sge_eqang(nstates,i)
          enddo
#ifdef _MPI_
          if ( iproc == 0 ) then
             write(*,'(a)') 'Bending reaction-coordinates in Multi-Window Tempering'
             write(*,1582)
             do i = 1, nbends_added_st
                idx = idx_bends(i)
                write(*,1583) atom_be1(idx), atom_be2(idx), atom_be3(idx), force_ang(idx), sge_eqang(1,i), sge_eqang(nstates,i)
             enddo
          endif
#endif
       endif

1582   format(5x, 'Atom1',5x,'Atom2',5x,'Atom3',5x,'K (kcal/mol/rad**2)',1x,'   alpha0 (degrees)',1x,'   alpha1 (degrees)')
1583   format( 3i10, 4x, 3(8x,f12.3) )

       if ( sge_tors ) then
          write(kprint,'(a)') 'Torsion reaction-coordinates in Multi-Window Tempering'
          write(kprint,1584)
          do i = 1, nitors_added_st
             idx = idx_tors(i)
             write(kprint,1585) atom_it1(idx), atom_it2(idx), atom_it3(idx), atom_it4(idx), &
                  force_died(idx), sge_eqdied(1,i), sge_eqdied(nstates,i)
          enddo
#ifdef _MPI_
          if ( iproc == 0 ) then
             write(*,'(a)') 'Torsion reaction-coordinates in Multi-Window Tempering'
             write(*,1584)
             do i = 1, nitors_added_st
                idx = idx_tors(i)
                write(*,1585) atom_it1(idx), atom_it2(idx), atom_it3(idx), atom_it4(idx), &
                     force_died(idx), sge_eqdied(1,i), sge_eqdied(nstates,i)
             enddo
          endif
#endif
       endif

1584   format(5x, 'Atom1', 5x, 'Atom2', 5x, 'Atom3', 5x, 'Atom4', 5x, &
               'K (kcal/mol/rad**2)', 1x, '    beta0 (degrees)', 1x, '    beta1 (degrees)' )
1585   format( 4i10, 4x, 3(8x,f12.3) )

    endif

    write(kprint,8007)
#ifdef _MPI_
    if ( iproc == 0 ) write(*,8007)
#endif

8100 FORMAT( // '=====  S E R I A L  G E N E R A L I Z E D  E N S E M B L E  S I M U L A T I O N  ====' / &
          /' ------ Number of replicas/processors:           ', i10 &
          /' ------ Number of ensembles:                     ', i10 &
          /' ------ Current replica/processor is:            ', i10 &
          /' ------ Ensemble of the current replica:         ', i10 &
          /' ------ Replica transition scheme:               ', 7x, a3 &
          /' ------ Replica transition attempts (frequency): ', f10.2, ' fs' &
          /' ------ Work updating (frequency):               ', f10.2, ' fs' &
          /' ------ PMF updating (frequency):                ', f10.2, ' fs' &
          /' ------ Number of samples in PMF update (0 means "all"):      ', i6 &
          /' ------ Print data for reweighting to the SGE_WHAM file (freq.): ', f10.2, ' fs' &
          /' ------ Print acceptance ratio to standard output (freq.):       ', f10.2, ' fs' )

#ifdef _MPI_
8101 FORMAT( // '====  S E R I A L  G E N E R A L I Z E D  E N S E M B L E  S I M U L A T I O N  ===='/ &
          /' ------ Number of replicas/processors:           ', i10 &
          /' ------ Number of ensembles:                     ', i10 &
          /' ------ Replica transition scheme:               ', 7x, a3 &
          /' ------ Replica transition attempts (frequency): ', f10.2, ' fs' &
          /' ------ Work updating (frequency):               ', f10.2, ' fs' &
          /' ------ PMF updating (frequency):                ', f10.2, ' fs' &
          /' ------ Number of samples in PMF update (0 means "all"):         ', i6 &
          /' ------ Print data for reweighting to the SGE_WHAM file (freq.): ', f10.2, ' fs' &
          /' ------ Print acceptance ratio to standard output (freq.):       ', f10.2, ' fs' )
#endif

8001 FORMAT( / ' ------ Local scaling/solute tempering is ON: definition of SOLUTE follows' )

8002 FORMAT( ' ------ group n. ',i4, ' goes from atom ', i6, ' to atom', i6 )

8003 FORMAT( ' ------ Total number of atoms of the solute is ',i6 )

8004 FORMAT( ' ------ Both solute-solvent and solute-solute are scaled' )

8005 FORMAT( ' ------ Only intrasolute interactions are scaled' )

8006 FORMAT( ' ------ Only solute-solvent interactions are scaled' )

8007 FORMAT( / '====================================================================================' / )

    return
  end subroutine sge_print_titles


!=============  INITIALIZE SGE_WHAM, SGE_ENERGY and SGE_DF FILES  ============
  subroutine init_files ( iproc, cpress )

    implicit none

    integer :: i
    integer, intent(in) :: iproc
    logical, intent(in) :: cpress

! initialization of SGE_WHAM file for WHAM application
    if ( sge_print_wham .and. (.not.always_accept)) then
       if ( sim_tempering_sge ) then
          write(ksge,187) &
               ( '      Bond index', '    Equil. value', '   Current value', i=1, nbonds_added_st ), &
               ( '      Bend index', '    Equil. value', '   Current value', i=1, nbends_added_st ), &
               ( '      Tors index', '    Equil. value', '   Current value', i=1, nitors_added_st )
187       format ( '#',3x,'Time(fs)',1x,'Ens.Index',1000(4x,a16) )
       else
          write(ksge,185)
185       format ( '#',3x,'Time(fs)',1x,'Ens.Index',5x,'Unscaled_Epot(1)',6x,'Unscaled_Epot(2)',6x,'Unscaled_Epot(3)' )
       endif
    endif

! initialization of SGE_ENERGY
    if ( sge_ts < 1.d-10 ) return

    open( unit=888, file='SGE_ENERGY', form='formatted' )
    if ( cpress ) THEN
       write(888,287)
287    format ( '#', 3x, 'Time(fs)', 1x, 'Ens.Index', 9x, 'Etot', 9x, 'Ekin', 9x, 'Epot', 8x, 'Ekinh', 8x, 'Epoth', &
            8x, 'EkinP', 11x, 'PV' )
    else
       write(888,288)
288    format ( '#', 3x, 'Time(fs)', 1x, 'Ens.index', 9x, 'Etot', 9x, 'Ekin', 9x, 'Epot', 8x, 'Ekinh', 8x, 'Epoth' )
    endif

! initialization of SGE_DF
    if(.not.always_accept) THEN 
    if ( sim_tempering .and. iproc == 0 ) then
       open( unit=889, file='SGE_DF', form='formatted' )
       write(889,'(a)') '# Df and Df_Error are dimensionless'
       write(889,289) '#   time(fs)', ( '      Df', i, ' -->', i+1, i=1, nstates-1), ( 'Df_Error', i, ' -->', i+1, i=1, nstates-1 )
289    format(a12,10000(2x,a8,i3,a4,i3))
    endif
    endif
          
    lfirst_sge = .false.

    return
  end subroutine init_files


!===============================  PRINT SGE_WHAM  ================================
  subroutine print_sge_wham ( time_fs, tors_added, eqdied, bend_added, eqang, bond_added, eqdist, radfact )

    use rem, only: para_index, rem_pot

    implicit none

    integer :: i
    real(8), intent(in) :: time_fs, radfact
    real(8), intent(in) :: eqdied(*), eqang(*), eqdist(*), tors_added(*), bend_added(*), bond_added(*)

    if ( sim_tempering_sge ) then
          write(ksge,'(f12.3,6x,i4,1x,1000(i20,2g20.8))') &
! time
               time_fs, &
! replica index
               para_index, &
! index, force constant, equilibrium and current value of BONDING reaction coordinates
               ( idx_bonds(i), eqdist(idx_bonds(i)), bond_added(idx_bonds(i)), i=1, nbonds_added_st ), &
! index, force constant, equilibrium and current value of BENDING reaction coordinates
               ( idx_bends(i), eqang(idx_bends(i)), bend_added(idx_bends(i))/radfact, i=1, nbends_added_st ), &
! index, force constant, equilibrium and current value of TORSION reaction coordinates
               ( idx_tors(i), eqdied(idx_tors(i)), tors_added(idx_tors(i))/radfact, i=1, nitors_added_st )
       else
          write(ksge,'(f12.3,6x,i4,1x,3(f20.8,2x))') &
! time
               time_fs, &
! replica index
               para_index, &
! unscaled potential
               rem_pot(1), rem_pot(2), rem_pot(3)
       endif

    return
  end subroutine print_sge_wham


!====================================  PRINT SGE_ENERGY  ===========================================
  subroutine print_sge_energy ( time_fs, ucek, pucek, ucepr, uceh, pext, volume, hpot, efact, un0, un1, unb, &
       rest_n1, rest_nb, cpress )

    use rem, only: para_index, rem_pot, rem_factor

    implicit none

    real(8) :: etot
    real(8), intent(in) :: time_fs, ucek, pucek, ucepr, uceh, pext, volume, hpot, efact, un0, un1, unb, rest_n1, rest_nb
    logical, intent(in) :: cpress

    etot = ( ucek + pucek + uceh + hpot * efact ) / 1000.d0

    if ( cpress ) etot = etot + ( ucepr + pext * volume * efact ) / 1000.d0

    if ( sim_tempering_sge ) then
       etot = etot + un0 + un1 + unb
    else
       etot = etot + rem_pot(1) * rem_factor(1) + rem_pot(2) * rem_factor(2) + rest_n1 + rem_pot(3) * rem_factor(3) + rest_nb
    endif

    if ( cpress ) then
       write(888,'(f12.3,6x,i4,1x,7(f12.3,1x))') &
! time
            time_fs, &
! replica index
            para_index, &
! total energy
            etot, &
! kinetic energy
            ( ucek + pucek ) / 1000.d0, &
! original (unscaled) potential energy
            un0 + un1 + unb, &
! thermostat kinetic energy
            uceh / 1000.d0, &
! thermostat potential energy
            hpot * efact / 1000.d0, &
! barostat kinetic energy
            ucepr / 1000.d0, &
! barostat potential energy
            ( pext * volume * efact ) / 1000.d0
    else
       write(888,'(f12.3,6x,i4,1x,5(f12.3,1x))') &
! time
            time_fs, &
! replica index
            para_index, &
! total energy
            etot, &
! kinetic energy
            ( ucek + pucek ) / 1000.d0, &
! original (unscaled) potential energy
            un0 + un1 + unb, &
! thermostat kinetic energy
            uceh / 1000.d0, &
! thermostat potential energy
            hpot * efact / 1000.d0
    endif

    return
  end subroutine print_sge_energy


!=========================    UPDATE DF   =============================
  subroutine sge_df_update ( time_fs, iproc )

    implicit none

! N_ATTEMPT: parameter which corresponds to the maximum number of failed attempts of free energy updates
! allowed before to come back to use the free-energy perturbation method (Jarzynski-like method)
! in place ofinstead of the BAR method (the value is arbitrary)
    integer, parameter :: n_attempt = 10

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(mpi_status_size), ierr
#endif

    integer :: i, j, nw_up, nw_down, nw_up_tot, nw_down_tot
    integer, intent(in) :: iproc
    real(8), intent(in) :: time_fs


! Update of free energy differences via Bennett (BAR) method
    do i = 1, nstates - 1

       nw_up = sge_up(i)
       nw_down = sge_down(i+1)

#ifdef _MPI_
       call mpi_allreduce ( nw_up, nw_up_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr )
       call mpi_allreduce ( nw_down, nw_down_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr )
#else
       nw_up_tot = nw_up
       nw_down_tot = nw_down
#endif
       if ( nw_down_tot >= nt_points_w .and. nw_up_tot >= nt_points_w ) then

          if ( estimate_df(i) == 0 .and. iproc == 0 ) then
             write(*,'(2a,i4,a,i4,a,f15.3,a)') 'Switch to Bennett '// &
                  'like method to update free energy difference',     &
                  ' Df ', i, ' -->', i+1, '  at time ', time_fs, ' fs'
          endif

          estimate_df(i) = 1
          estimate_df_up(i) = 0
          estimate_df_down(i) = 0

! stima DF con metodo di Bennett
          call bnt ( nw_up_tot, nw_down_tot, nw_up, nw_down, w_up(1,i), w_down(1,i+1), dfi(i), df_var(i) )

          df_av(i)  = df_av(i)  + dfi(i) / df_var(i)
          ndf_av(i) = ndf_av(i) + 1.d0   / df_var(i)

          if ( iproc == 0 .and. sge_update == 0 .and. sge_update_fly <= 0 ) then
             write(*,'(2(a,i4))') 'All samples are used to calculate '// &
                  'weighted free energies. Updating Df', i, ' -->', i+1
          endif

! Update using a limited number of Df estimates (from input)
          if ( sge_update > 0 ) then
             if ( iproc == 0 .and. sge_update_fly <= 0 )  write(*,'(a,1x,i10,3x,2(a,i4))')      &
                  'Number of samples used to compute weighted free '//'energies (from input):', &
                  sge_update, 'Updating Df', i, ' -->', i+1
             call sge_df_update_1 (i)
          endif

! Update on-the-fly by using a limited number of Df estimates (from files SGE_DF_FLY.set and SGE_DF_FLY.dat)
          if ( sge_update_fly > 0 ) then
             if ( iproc == 0 ) then
                write(*,'(a,1x,i10,3x,2(a,i4))') 'Number of samples '//                                          &
                     'used to compute weighted free energies '//'(change on-the-fly from file SGE_DF_FLY.dat):', &
                     sge_update_fly, 'Updating Df', i, ' -->', i+1
                nfly_df = nfly_df + 1
                write(444,*) i, dfi(i), df_var(i)
                call sge_df_update_2 (i)
             endif
#ifdef _MPI_
             call mpi_bcast( df(i), 1, mpi_double_precision, 0, mpi_comm_world, ierr )
             call mpi_bcast( ndf_av_fly(i), 1, mpi_double_precision, 0, mpi_comm_world, ierr )
#endif
! Update by using either a limited number of Df estimates (from input: sge_update > 0)
! or all Df estimates (sge_update = 0)
          else
             if ( ndf_av(i) /= 0.d0 ) df(i) = df_av(i) / ndf_av(i)
          endif

! Work counters are zeroed
          sge_up(i) = 0
          sge_down(i+1) = 0

       elseif ( nw_up_tot < nt_points_w .and. nw_down_tot >= 2*nt_points_w ) then ! the downward work matrix is reduced in size
                                                                                  ! in order to keep the numbers of upward and downward
                                                                                  ! work samples comparable each to other
          sge_down(i+1) = nint( sge_down(i+1) / 2. )
          do j = 1, sge_down(i+1)
             w_down(j,i+1) = w_down(j+sge_down(i+1),i+1)
          enddo

       elseif ( nw_down_tot < nt_points_w .and. nw_up_tot >= 2*nt_points_w ) then ! the upward work matrix is reduced in size
                                                                                  ! in order to keep the numbers of upward and downward
                                                                                  ! work samples comparable each to other
          sge_up(i) = nint( sge_up(i) / 2. )
          do j = 1, sge_up(i)
             w_up(j,i) = w_up(j+sge_up(i),i)
          enddo

       endif

! Return back to Jarzynski method (from Bennett method) for updating free
! energies if after n_attempt attempts the free energy is not updated
       if ( estimate_df(i) == 1 ) then

          if ( nw_down_tot == nw_down_tot_old(i) ) then
             nconst_down(i) = nconst_down(i) + 1
          else
             nconst_down(i) = 0
             nw_down_tot_old(i) = nw_down_tot
          endif

          if ( nconst_down(i) == n_attempt ) then
             if ( iproc == 0 ) write(*,'(2a,i4,a,i4,a,f15.3,a)') 'Switch back to free-energy-perturbation '// &
                  'method to update free energy difference', ' Df ', i, ' -->', i+1, '  at time ', time_fs, ' fs'
             nconst_down(i) = 0
             nw_down_tot_old(i) = nw_down_tot
             estimate_df(i) = 0
             estimate_df_up(i) = 1
             estimate_df_down(i) = 1
! The weighted averages of the free energy differences can be zeroed
!             ndf_av(i) = 0.
!             df_av(i) = 0.
! Work counters are zeroed
             sge_up(i)= 0
             sge_down(i+1) = 0
          endif

          if ( nw_up_tot == nw_up_tot_old(i) ) then
             nconst_up(i) = nconst_up(i) + 1
          else
             nconst_up(i) = 0
             nw_up_tot_old(i) = nw_up_tot
          endif

          if ( nconst_up(i) == n_attempt ) then
             if ( iproc == 0 ) write(*,'(2a,i4,a,i4,a,f15.3,a)') 'Switch back to free-energy-perturbation '// &
                  'method to update free energy difference', ' Df ', i, ' -->', i+1, '  at time ', time_fs, ' fs'
             nconst_up(i) = 0
             nw_up_tot_old(i) = nw_up_tot
             estimate_df(i) = 0
             estimate_df_up(i) = 1
             estimate_df_down(i) = 1
! The weighted averages of the free energy differences can be zeroed
!             ndf_av(i) = 0.
!             df_av(i) = 0.
! Work counters are zeroed
             sge_up(i) = 0
             sge_down(i+1) = 0
          endif

       endif



! Update of free energy differences via Jarzynski method
       if ( estimate_df(i) == 0 ) then

          nw_up = sge_up_jar(i)
          nw_down = sge_down_jar(i+1)

#ifdef _MPI_
          call mpi_allreduce ( nw_up, nw_up_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr )
          call mpi_allreduce ( nw_down, nw_down_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr )
#else
          nw_up_tot = nw_up
          nw_down_tot = nw_down
#endif
! DF = F(T_i+1) - F(T_i) = ln{ Z(T_i) / Z(T_i+1) }        (dimensionless)
          if ( nw_down_tot >= nt_points_w ) then
             estimate_df_down(i) = 1
             df_down(i) = log( jarzy(w_down_jar(1,i+1),nw_down) / nw_down_tot )
! Work counters are zeroed
             sge_down_jar(i+1) = 0
          elseif ( nw_up_tot >= nt_points_w ) then
             estimate_df_up(i) = 1
             df_up(i) = -log( jarzy( w_up_jar(1,i), nw_up ) / nw_up_tot )
! Work counters are zeroed
             sge_up_jar(i) = 0
          endif

       endif


! Check if the variables ESTIMATE_DF* are consistent
       if ( estimate_df(i) == 1 .and. ( estimate_df_up(i) == 1 .or. estimate_df_down(i) == 1 ) ) then
          write(*,'(a)') 'Variables ESTIMATE_DF* are not consistent.'
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort( mpi_comm_world, ierr )
#else
          stop
#endif
       endif

    enddo

    return
  end subroutine sge_df_update



!================   Update free energy with last Df estimates  ============
  subroutine sge_df_update_1 ( i )

    implicit none

    integer, intent(in) :: i

    cnt_upd(i) = cnt_upd(i) + 1
    if ( cnt_upd(i) > sge_update ) cnt_upd(i) = 1
    if ( df_var_inst(i,cnt_upd(i)) /= 0.d0 ) then
       df_av(i) = df_av(i) - dfi_inst(i,cnt_upd(i)) / df_var_inst(i,cnt_upd(i))
       ndf_av(i) = ndf_av(i) - 1.d0 / df_var_inst(i,cnt_upd(i))
    endif
    dfi_inst(i,cnt_upd(i)) = dfi(i)
    df_var_inst(i,cnt_upd(i)) = df_var(i)

    return
  end subroutine sge_df_update_1

!======   Update free energy with last Df estimates (on the fly)  ========
  subroutine sge_df_update_2 ( i )

    implicit none

    integer :: ns, nfly_flag, nfly_count, k, i
    real(8) :: df_av_tmp, ndf_av_tmp, dfi_tmp, df_var_tmp

    rewind(444)
    nfly_count = 0
    do k = 1, nfly_df
       read(444,*) ns
       if ( ns == i ) nfly_count = nfly_count + 1
    enddo

    rewind(444)
    nfly_flag = nfly_count - sge_update_fly
    df_av_tmp = 0.d0
    ndf_av_tmp = 0.d0
    nfly_count = 0

    do k = 1, nfly_df
       read(444,*) ns, dfi_tmp, df_var_tmp
       if ( ns == i ) then
          nfly_count = nfly_count + 1
          if ( nfly_count > nfly_flag ) then
             df_av_tmp  = df_av_tmp  + dfi_tmp / df_var_tmp
             ndf_av_tmp = ndf_av_tmp + 1.d0    / df_var_tmp
          endif
       endif
    enddo

    if ( ndf_av_tmp /= 0.d0 ) then
       df(i) = df_av_tmp / ndf_av_tmp
       ndf_av_fly(i) = ndf_av_tmp
    endif

    return
  end subroutine sge_df_update_2





!=============  UPDATE SIZE OF WORK ARRAYS  =================
  subroutine sge_size_update ( iproc )

    use rem, only: para_index

    implicit none

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(mpi_status_size), ierr
#endif

    real(8), allocatable :: wtmp(:,:)
    integer :: dim, i
    integer, intent(in) :: iproc

! Check and possible increase of the matrix W_UP()
    dim = sge_up(para_index) + 1

    if ( dim > current_dim_up ) then
       allocate( wtmp(current_dim_up,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for WTMP()'
       wtmp(1:current_dim_up,:) = w_up(1:current_dim_up,:)

       deallocate( w_up, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for W_UP()'

       allocate( w_up(dim,nstates), STAT=IERROR )
       if ( ierror /= 0 ) stop 'Allocation failed for W_UP()'
       w_up(:,:) = 0.d0
       w_up(1:current_dim_up,:) = wtmp(1:current_dim_up,:)

       deallocate( wtmp, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for WTMP()'

       current_dim_up = dim
    endif

! Check and possible increase of the matrix W_DOWN()
    dim = sge_down(para_index) + 1

    if ( dim > current_dim_down ) then
       allocate( wtmp(current_dim_down,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for WTMP()'
       wtmp(1:current_dim_down,:) = w_down(1:current_dim_down,:)

       deallocate( w_down, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for W_DOWN()'

       allocate( w_down(dim,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for W_DOWN()'
       w_down(:,:) = 0.d0
       w_down(1:current_dim_down,:) = wtmp(1:current_dim_down,:)

       deallocate( wtmp, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for WTMP()'

       current_dim_down = dim
    endif

! Check dimensions
    if ( para_index /= nstates ) then
       if ( sge_up(para_index) + 1 > current_dim_up ) then
          write(*,'(a)') 'W_UP is undersized'
          write(*,'(a,1x,i8)') 'Current dimension:', sge_up(para_index)
          write(*,'(a,1x,i8)') 'Maximum dimension:', current_dim_up
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
    endif

    if ( para_index /= 1 ) then
       if ( sge_down(para_index) + 1 > current_dim_down ) then
          write(*,'(a)') 'W_DOWN is undersized'
          write(*,'(a,1x,i8)') 'Current dimension:', sge_down(para_index)
          write(*,'(a,1x,i8)') 'Maximum dimension:', current_dim_down
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
    endif

! Check and possible increase of the matrix W_UP_JAR()
    dim = sge_up_jar(para_index) + 1
              
    if ( dim > current_dim_up_jar ) then
       allocate( wtmp(current_dim_up_jar,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for WTMP()'
       wtmp(1:current_dim_up_jar,:) = w_up_jar(1:current_dim_up_jar,:)

       deallocate( w_up_jar, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for W_UP_JAR()'

       allocate( w_up_jar(dim,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for W_UP_JAR()'
       w_up_jar(:,:) = 0.d0
       w_up_jar(1:current_dim_up_jar,:) = wtmp(1:current_dim_up_jar,:)

       deallocate( wtmp, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for WTMP()'

       current_dim_up_jar = dim
    endif

! Check and possible increase of the matrix W_DOWN_JAR()
    dim = sge_down_jar(para_index) + 1

    if ( dim > current_dim_down_jar ) then
       allocate( wtmp(current_dim_down_jar,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for WTMP()'
       wtmp(1:current_dim_down_jar,:) = w_down_jar(1:current_dim_down_jar,:)

       deallocate( w_down_jar, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for W_DOWN_JAR()'

       allocate( w_down_jar(dim,nstates), stat=ierror )
       if ( ierror /= 0 ) stop 'Allocation failed for W_DOWN_JAR()'
       w_down_jar(:,:) = 0.d0
       w_down_jar(1:current_dim_down_jar,:) = wtmp(1:current_dim_down_jar,:)

       deallocate( wtmp, stat=ierror )
       if ( ierror /= 0 ) stop 'Deallocation failed for WTMP()'

       current_dim_down_jar = dim
    endif

! Check dimensions
    if ( para_index /= nstates ) then
       if ( sge_up_jar(para_index) + 1 > current_dim_up_jar ) then
          write(*,'(a)') 'W_UP_JAR is undersized'
          write(*,'(a,1x,i8)') 'Current dimension:', sge_up_jar(para_index)
          write(*,'(a,1x,i8)') 'Maximum dimension:', current_dim_up_jar
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
    endif

    if ( para_index /= 0 ) then
       if ( sge_down_jar(para_index) + 1 > current_dim_down_jar ) then
          write(*,'(a)') 'W_DOWN_JAR is undersized'
          write(*,'(a,1x,i8)') 'Current dimension:', sge_down_jar(para_index)
          write(*,'(a,1x,i8)') 'Maximum dimension:', current_dim_down_jar
          write(*,'(a)') 'The program will stop.'
#ifdef _MPI_
          call mpi_abort ( mpi_comm_world, ierr )
#endif
          stop
       endif
    endif

    return
  end subroutine sge_size_update








!==================  FUNCTIONS AND SUBROUTINES USED IN CALCULATING DF  ==================

! Compute the free energy difference DA using the Bennett method
  subroutine bnt ( nwf_tot, nwb_tot, nwf, nwb, wf, wb, da, var_da )

    implicit none

    integer, intent(in) :: nwf, nwb, nwf_tot, nwb_tot
    real(kr8), intent(in) :: wf(nwf), wb(nwb)
    real(kr8), intent(out) :: da, var_da
    real(kr8) :: da1

! The adimensional free energy difference is calculated
    da1 = rtsafe(nwf,nwb,wf,wb,-50000.0_kr8,50000.0_kr8,0.0001_kr8)
    da = da1 - log( dble(nwb_tot) / dble(nwf_tot) )

! The variance of the adimensional free energy difference is calculated using Bennett method (Phys. Rev. Lett. n140601 v91 y2003)
    var_da = variance_bnt ( nwf_tot, nwb_tot, nwf, nwb, wf, wb, da1 )

    return
  end subroutine bnt


!--------------------------------------------------------------------
  function rtsafe ( nwf, nwb, wf, wb, x1, x2, xacc )

    implicit none

    integer, parameter :: maxit = 100
    integer, intent(in) :: nwf, nwb
    real(kr8), intent(in) :: wf(nwf), wb(nwb)
    real(kr8) :: rtsafe, x1, x2, xacc
    integer :: j
    real(kr8) :: df, dx, dxold, f, fh, fl, temp, xh, xl
    
    call funcd( nwf, nwb, wf, wb, x1, fl, df )
    call funcd( nwf, nwb, wf, wb, x2, fh, df )

    if ( ( fl > 0.0_kr8 .and. fh > 0.0_kr8 ) .or. ( fl < 0.0_kr8 .and. fh < 0.0_kr8 ) ) &
         write(*,*) 'Root must be bracketed in rtsafe'
    if ( fl == 0.0_kr8 ) then
       rtsafe = x1
       return
    elseif ( fh == 0.0_kr8 ) then
       rtsafe = x2
       return
    elseif ( fl < 0.0_kr8 ) then
       xl = x1
       xh = x2
    else
       xh = x1
       xl = x2
    endif
    rtsafe = 0.5_kr8 * ( x1 + x2 )
    dxold = abs( x2 - x1 )
    dx = dxold
    call funcd( nwf, nwb, wf, wb, rtsafe, f, df )
    do j = 1, maxit
       if ( ( (rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.0_kr8 .or. abs(2.*f) > abs(dxold*df) ) then 
          dxold = dx
          dx = 0.5_kr8 * ( xh - xl )
          rtsafe = xl + dx
          if ( xl == rtsafe ) return
       else
          dxold = dx
          dx = f / df
          temp = rtsafe
          rtsafe = rtsafe - dx
          if ( temp == rtsafe ) return
       endif
       if ( abs(dx) < xacc ) return
       call funcd( nwf, nwb, wf, wb, rtsafe, f, df )
       if ( f < 0.0_kr8 ) then
          xl = rtsafe
       else
          xh = rtsafe
       endif
    enddo
    write(*,*) 'rtsafe exceeding maximum iterations'

    return
  end function rtsafe

!--------------------------------------------------------------------
  subroutine funcd ( nwf, nwb, wf, wb, da, ll, lld )

    implicit none

    real(kr8), parameter :: eps = 0.00001_kr8
    integer, intent(in) :: nwf, nwb
    real(kr8), intent(in) :: wf(nwf), wb(nwb)
    real(kr8), intent(in) :: da
    real(kr8), intent(out) :: ll, lld

    ll = dllike( nwf, nwb, wf, wb, da )
    lld = ( dllike( nwf, nwb, wf, wb, da + eps ) - dllike( nwf, nwb, wf, wb, da - eps ) ) / ( 2.0_kr8 * eps )
    
  end subroutine funcd

!--------------------------------------------------------------------
  function dllike ( nwf, nwb, wf, wb, da )

    implicit none

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(mpi_status_size), ierr
    real(kr8) :: ft
#endif
    integer, intent(in) :: nwf,nwb
    real(kr8), intent(in) :: wf(nwf), wb(nwb)
    real(kr8), intent(in) :: da
    real(kr8) :: f
    real(kr8) :: dllike

    f = sum( 1.0_kr8 / (1.0_kr8 + exp(wf - da)) ) - sum( 1.0_kr8 / (1.0_kr8 + exp(wb + da)) )
#ifdef _MPI_
    call mpi_allreduce ( f, ft, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
    dllike = ft
#else
    dllike = f
#endif

  end function dllike

!--------------------------------------------------------------------
! The variance VAR of the free energy difference DA is computed
  function variance_bnt ( nwf_tot, nwb_tot, nwf, nwb, wf, wb, da1 )

    implicit none

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(mpi_status_size), ierr
    real(8) :: ft
#endif
    integer, intent(in) :: nwf, nwb, nwf_tot, nwb_tot
    real(8), intent(in) :: wf(nwf), wb(nwb), da1
    real(8) :: f, variance_bnt

    f = sum( 1. / ( 1. +  cosh(wf(:) - da1) ) ) + sum( 1. / ( 1. + cosh(WB(:) + DA1) ) )
#ifdef _MPI_
    call mpi_allreduce ( f, ft, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
    variance_bnt = 2. / ft - 1. / dble(nwf_tot) - 1. / dble(nwb_tot)
#else
    variance_bnt = 2. / f  - 1. / dble(nwf_tot) - 1. / dble(nwb_tot)
#endif
    
  end function variance_bnt

!--------------------------------------------------------------------
! Free energy difference is calculated using Jarzynski-like method
  function jarzy ( w, nw )

    implicit none

#ifdef _MPI_
    include 'mpif.h'
    integer :: status(mpi_status_size), ierr
    real(8) :: ft
#endif
    integer :: nw
    real(8) :: w(nw), jarzy, f

    f = sum( exp(-w(:)) )
#ifdef _MPI_
    call mpi_allreduce ( f, ft, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
    jarzy = ft
#else
    jarzy = f
#endif

    return
  end function jarzy

end module sge
