SUBROUTINE blkdta
!************************************************************************
!*                                                                      *
!*     Subroutine BLKDTA set up defaults for the simulation             *
!*     parameters.                                                      *
!*                                                                      *
!************************************************************************

  use unit
  use parst
  use cpropar
  use meta
  use rem
  use sge, only: sge_update, sge_update_rst, sge_update_fly,nstep_fly,lfirst_sge,sge_tsave_work, &
       sge_tsave_pmf,sim_tempering, sge_ts, sim_tempering_sge, sge_print_wham,sge_print, sge_f,  & 
       sge_transition_scheme,always_accept,ncycle,ncycle_reject,use_saved,go_save,nfreqsge, & 
       print_always
  use dcd, only: kdcd,ndcd,dcd_frame,dcd_noh
#ifdef _MPI_
  use orac_mpi
#endif
  use omp_integr   ! this is for setting up times in NON-OpenMP version

  IMPLICIT none
  
  INTEGER i,j
  
  COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
  INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)
  CHARACTER*80 fileout

!=======================================================================
!*------- Set defaults parameters ---------------------------------------
!*=======================================================================

#ifdef _MPI_
  WRITE (fileout,'(a,i4.4)') TRIM("out"),iproc_t ! sets the PARXX dirs to iproc_t.  
  kprint=85
  open(unit=kprint,file=fileout,form='FORMATTED',status='UNKNOWN')
  !       this variable controls who is entitled to print in the traj
  !       communication group in normal execution.
  knlist=86
  open(unit=knlist,file='INPUT',form='FORMATTED',status='UNKNOWN')
#else
  knlist=5
  kprint=6
#endif
! with SAVE_ALL_RESTARTS, kdump is not incremented by openf.f
! first used value kdump=0 will result in an error in AIX 6.1 (UNIT=0 is stderr)
!      kdump=0
  ffield_type="amber" ! this is the default force field
  nthr =1
  nthr1=1
  nthr2=1
  indtar=.false.
  indlig=.false.
  strbonds=0
  strbends=0
  strtors=0
  time_n0=0.d0
  time_n1=0.d0
  time_h=0.d0
  time_l=0.d0
  time_m=0.d0
  time_u=0.d0
  time_f=0.d0
  time_cn=0.d0
  time_in=0.d0
  tcpu_n0=0.d0
  tcpu_n1=0.d0
  tcpu_h=0.d0
  tcpu_l=0.d0
  tcpu_m=0.d0
  tcpu_u=0.d0
  tcpu_f=0.d0
  tcpu_cn=0.d0
  tcpu_in=0.d0
  kdump=10
  kgdata=0
  kcnfig=0
  kgofr=0
  kplot=0
  kconf=0
  print_qwork=.false.
  itprint=1
  ifprint=1
  kpot=0
  kout=0
  kfield=0
  kvaf=0
  krms=0
  kxrms=0
  kxrms_atm=0
  ksol1=0
  ksol2=0
  kgroup=0
  ktop=0
  ktemplate=0
  ktest=0
  iurest=0
  kspec=0
  kgr=0
  kplot_fragm=0
  kplot_center=0
  kplot_steer = 0
  kplot_stn0=0
  kgofr_sk=0
  kcoord_slt=0
  kcoord_slv=0
  ktpgprm_read=0
  ktpgprm_write=0
  kfreq=0
  annpro=0
  jerr=1
  aaxis=500.0
  baxis=500.0
  caxis=500.0
  alf=90.0
  bet=90.0
  gam=90.0
  frequencies=.FALSE.
  debug=.FALSE.
  debug_rs=.FALSE.
  debug_ct=.FALSE.
  debug_st=.FALSE.
  debug_bt=.FALSE.
  debug_pt=.FALSE.
  debug_it=.FALSE.
  check_coord=.TRUE.
  do_workq=.false.
  energy_then_die=.FALSE.
  docking=.FALSE.
  cutoff=.TRUE.
  grpcut=.FALSE.
  protei=.FALSE.
  solven=.FALSE.
  electr=.FALSE.
  clewld=.FALSE.
  linser=.FALSE.
  recstrc=.FALSE.
  stinit = .FALSE.
  tpgfil = .FALSE.
  rescm  = .FALSE.
  tpgwbn = .FALSE.
  prttpg = .FALSE.
  sction=.FALSE.
  prtatl=.FALSE.
  prthbl=.FALSE.
  prtnbl=.FALSE.
  prtn14=.FALSE.
  prtcnl=.FALSE.
  prtbal=.FALSE.
  prtptl=.FALSE.
  prtitl=.FALSE.
  prtipw=.FALSE.
  hmass=.FALSE.
  prttor=.FALSE.
  prtben=.FALSE.
  prtstr=.FALSE.
  path_steer=.FALSE.
  igr=.FALSE.
  gofr=.FALSE.
  restart_cont=.FALSE.
  read_co=.FALSE.
  gofr_nprint=0
  gofr_navg=0
  gofr_ncomp=0
  gofr_neighbor=.TRUE.
  gofr_intra=.FALSE.
  gofr_avg=.FALSE.
  gofr_cut=0.0D0
  delrg=0.05D0
  lenerg=.FALSE.
  uenerg=.FALSE.
  fold=.FALSE.
  wrtgyr=.FALSE.
  noneq=.FALSE.
  lfield=.FALSE.
  print_unformatted=.FALSE.
  scharge=.FALSE.
  nonbnd=.TRUE.
  hoover=.FALSE.
  thermos=.FALSE.
  vacf=.FALSE.
  mdsim=.TRUE.
  minimize=.FALSE.
  write_grad=.FALSE.
  annealing=.FALSE.
  change_cell=.FALSE.
  analys=.FALSE.
  anxrms=.FALSE.
  anxrms_cell=.FALSE.
  ascii_nocell=.FALSE.
  xyz_nocell=.FALSE.
  anxca=.FALSE.
  anxbc=.FALSE.
  anxhe=.FALSE.
  anrms=.FALSE.
  anslv=.FALSE.
  stretch=.FALSE.
  stretch_heavy=.FALSE.
  stoprun=.FALSE.
  md_respa=.FALSE.
  adihed=.TRUE.
  hydbnd=.FALSE.
  spectra=.FALSE.
  lgr=.FALSE.
  linked_cell=.FALSE.
  agbnp=.false.
  lsasa=.false.
  restart_cont=.FALSE.
  adjust_cnstr=.TRUE.
  slv_read=.FALSE.
  slv_exist=.FALSE.
  slt_exist=.FALSE.
  slv_generate=.FALSE.
  slv_create=.FALSE.
  pdb_read=.FALSE.
  slt_create=.FALSE.
  abmd=.FALSE.
  dissociate=.FALSE.
  associate=.FALSE.
  diss_mol=.FALSE.
  diss_atoms=.FALSE.
  abmd_tors=.FALSE.
  occupy_space=.FALSE.
  diffusion=.FALSE.
  time_corr=.FALSE.
  not_time_corr=.FALSE.
  erf_corr=.FALSE.
  erfc_spline=.FALSE.
  remove_momentum=.FALSE.
  always_scale_nose=.FALSE.
  anprot=.FALSE.
  addstrcom=.FALSE.
  steer_com=.FALSE.
  bussi=.FALSE.
  moveitors=.TRUE.
  gen_itp=.FALSE.
  read_remark=.false.
  divide_records=1
  divide_spline=1
  ns1=256
  nv1=256
  ns2=512
  nv2=512
  ns3=1024
  nv3=1024
  ns4=2048
  nv4=2048
  niop=1
  nios=1
  rpass =0.005
  ngprint = 10000
  nprot_charges=0
  corr_atoms(1)=0
  do i=1,5
     ntmtsp(i)=-1
  end do
  do i=1,3
     ntmtss(i)=-1
  end do

  totvaf=300
  novaf=3
  pext=0.0D0
  taut=0.0D0
  taup=0.0D0
  compressibility=5.3D-04
  isostress=.FALSE.
  isostressxy=.FALSE.
  cpress=.FALSE.
  pressure=.FALSE.
  lberendsen=.FALSE.
  landersen=.FALSE.
  coupl_mol=.TRUE.
  coupl_grp=.FALSE.
  coupl_atm=.FALSE.
  wpr=1.0D19
  massinfty=1.0D20
  ecut=1.0D0
  epcut=1.0D0
  echrg=-11.0D0
  pec=8
  nec=100
  rspcut=1.0D0
  rspoff=0.0D0
  rspon=1.0D0
  hrcut=1.0D0
  hacut=10.0D0
  hrson=3.50D0
  hrsoff=4.50D0
  hanon=50.0D0
  hanoff=70.0D0
  rcut_hb=4.5D0
  acut_hb=0.0D0
  a2cut_hb=0.0D0
  dock_radius=0.d0
  vsite=1660.d0
  nhskip=3
  nupdte=20
  nupdte_index=1
  alphal=0.4D0
  omega=1.0D0
  spring=0.0D0
  lmbdt=0.5D0
  rgnato=0
  fudge=1.0D0
  lj_fudge=1.0D0
  dock_radius=0.d0
  ancorr=100
  antstar=0.0D0
  anfact=1.0D0
  dtemph=5000.0D0
  dtemppr=10000.0d0 
  hstep_freq=0.03D0
  nstep_freq=6
  neta=3
  DO i=1,neta
     qmass(i)=0.0D0
  END DO
  scale=1E8
  iseed=1
  compressibility=4.9D-04
  eps_energy=0.5D0
  itor_ptype=0
  nospt=3
  nplot_fragm=0
  nplot_center=0
  ncofactor=0
  nfragm=0
  nrespa=1
  n0respa=1
  n1respa=1
  mrespa=1
  lrespa=1
  ltest_times=.false.
  clean=.true.
  start_conf= .false.
  distmax = 0.25
  rcuth=9.7
  rtolh=0.3
  rcutl=7.3
  rtoll=0.3
  rcutm=4.1
  rtolm=0.3
  rneih=1.5 
  rneil=0.45
  rneim=0.35
  kl=0.0
  klt=0.0
  km=0.0
  kmt=0.0
  erfc_bin=0.007D0
  pme=.false.
  shell_pme='0'
  nfft1=0
  nfft2=0
  nfft3=0
  ntimes=0
  mim_lim=20
  pme_order=0
  DO  i=1,8
     DO  j=1,3
        rmol(j,i)=0.0D0
     END DO
  end do
  mixrww=.true.
  mixrwp=.false.
  amphi=.FALSE.
  bending=.TRUE.
  replicate=.FALSE.
  write_ff_pars=.FALSE.
  ndocking=0
  jligand=-1
    nplot=0
  nascii=0
  nfield=0
  mfield=1
  nvi = 0
  ndipole = 0
  nhoov=1
  whoov=500.0D0
  time=1.0D0
  radius=1.6d+0
  rkcut=1.0D9
  rcut_corr = 0.0d0 
  nflag(1)=0
  nflag(2)=0
  maxstp=1
  maxrun=0
  nrject=1
  nprop=200
  nsave=0
  nconf=0
  nprint=1
  menerg=0
  t=300.0D0
  nform=1
  icl=1
  icm=1
  icn=1
  nmol=0
  nato_slv=0
  start_anl=1
  stop_anl=2**(24)
  update_anl=0
  hbonds_res=.FALSE.
  hbonds_tot=.FALSE.
  prttopl=.FALSE.
  fragm_dist=.FALSE.
  ! Statements for REM
#ifdef _MPI_
  use_saved=.false.
  go_save=.false.
#endif
  l_mpi=.false.
  rem_run = .false.
  ncycle=0
  rem_acc=0
  rep_acc=0
  rem_ts=0.d0
  rem_print=0.d0
  rem_segment = .false.
  prnt_remdt = .false.
  lskip_force= .false. 
  nrem = 0
  nbatteries=1
  always_accept=.false.
  ! Statements for SGE
  print_always=0
  nfreqsge=0
  sge_ts = 0.d0
  sge_f = 0.d0
  sge_update = 0
  sge_update_rst = 0
  sge_update_fly = 0
  nstep_fly = 0
  sge_tsave_work = 0.d0
  sge_tsave_pmf = 0.d0
  sge_print = 0.d0
  sim_tempering = .false.
  sim_tempering_sge = .false.
  lfirst_sge = .true.
  sge_print_wham = .false.
  sge_transition_scheme = 'DEO'
  ! Statements for metadynamics
  kmeta = 0
  meta_run = .FALSE.
  meta_tempered = .FALSE.
  meta_wtempered = .FALSE.
  metar = .FALSE.
  metaw = .FALSE.
  meta_analysis = .FALSE.
  meta_ts = 0.d0
  meta_height = 0.d0
  meta_nrc = 0
  meta_print = 0.d0
  meta_n = 0
  meta_step = 0
  meta_tot = 0.d0
  meta_maxv = 0.d0
  meta_temperature = 0.d0
  ! statements for dcd
  kdcd = 0
  krem_dcd = 0
  ndcd = 0
  dcd_frame = 0
  dcd_noh = .false.
  omp_timing=.false.
  ! statements for alchemy
  alchemy=.FALSE.
  intramol=.false.
  ljintra=.false.
  shifted=.false.
  wrk_alchemy=0.d0
  qsoft=0.5   ! Beutler sof-core gamma parameter for Coulomb
  ljsoft=0.5  ! Beutler sof-core gamma parameter for LJ
  qshift=4.0
  ljshift=0.35
#ifdef _OMP_
  omp_dynamic=.false.
  level_integr=1
  level_cnst=1
  level_n0=1
  level_n1=1
  level_shell_m=1
  level_shell_l=1
  level_update=1
  level_fftw=1
  cache_line_size=8
  cache_line_size1=8
  cache_line_size2=8
! initialize logical array for _OMP_ to avoid false sharing 
  array_omp=.false.
  array1_omp=.false.
  array2_omp=.false.
  friction=0.1
#endif

!================= END OF EXECUTABLE STATEMENTS ========================

RETURN
END SUBROUTINE blkdta

