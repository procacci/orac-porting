module meta
  implicit none
  ! metadim_nhill : max number of hills 
  ! metadim_nrc   : max number of reaction coordinates
  ! metadim_nat   : max number of atoms involved in the definition of a single reaction coordinate
  ! temporarily reduced to save memory allocation
  !  integer, parameter :: metadim_nhill=1000000, metadim_nrc=10, metadim_nat=100
  integer, parameter :: metadim_nhill=20000, metadim_nrc=2, metadim_nat=100
  integer :: meta_n                               ! current number of hills
  integer :: meta_step                            ! current number of steps
  integer :: meta_nrc                             ! total number of reaction coordinates (RCs)
  integer :: kmeta,kmetar,kmetaw                  ! output units
  integer :: meta_atoms(metadim_nat,metadim_nrc) ! each row contains the atoms that define a RC
  real(8) :: meta_ts           ! hills deposition interval 
  real(8) :: meta_print                           ! printing interval
  real(8) :: meta_height                          ! hill heigth
  real(8) :: meta_tot          ! sum_i 1/V_i over the hills
  real(8) :: meta_maxv         ! max value of the biasing potential at current time
  real(8) :: meta_temperature  ! temperature for the tempered-well tempered metadynamics algorithm
  real(8) :: meta_rc(metadim_nrc)                 ! current position of the system in the RCs space
  real(8) :: meta_width(metadim_nrc)              ! array of the widths of an hill along each RC
  real(8) :: meta_e(metadim_nhill)                ! array of the biasing potential corresponding to the center of each of the hills
  real(8) :: meta_list(metadim_nhill,metadim_nrc) ! the i-th row is the position of the i-th hill in the RCs space
  real(8) :: meta_weight(metadim_nhill) ! array of weights for well-tempered metadynamics
  logical :: meta_run                                 ! if true, the program performs a metadynamics run
  logical :: meta_analysis     ! if true, the program performs a metadynamics run
  logical :: meta_tempered   ! if true, tempered metadynamics
  logical :: meta_wtempered   ! if true, well tempered metadynamics
  logical :: metar,metaw                          ! if true, read/print a restart file
  character(len=1)  :: meta_dc                          ! dummy character
  character(len=4)  :: meta_rctype(metadim_nrc)         ! array containing the kinds of the RCs  
  character(len=80) :: metar_file,metaw_file,meta_data_file
  
end module meta

