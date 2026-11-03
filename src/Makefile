#
# Makefile for ORAC generated with configure script  
#    Target : GNU-FFTW-OMP-MPI
# =========== non specific variables ==========================


.DEFAULT_GOAL := GNU-FFTW-OMP-MPI
# look for prerequisites in parent dir and other dirs:
VPATH =..:../libnblist:../libagbnp:/usr/local/lib
#
CPP = /usr/bin/cpp -traditional
CFLAGS = -O3 -I../libnblist -I../libagbnp
# uncomment this if you want printout in the agbnp routine
#CFLAGS := -DAGBNP_WRITE_APARAMS -DAGBNP_WRITE_STRUCT $(CFLAGS)

# if fortran compiler underscores symbols and default is not working
#   (linker gives an error like "undefined reference to `vstart_'")
# then uncomment next line and/or add "-DFC_UNDERSCORES" to 
#   an explicitly set CFLAGS in your target
# CFLAGS := $(CFLAGS) -DFC_UNDERSCORES

FFLAGS =  -I.. -I. -I/usr/local/include -I/usr/include -O3
LDFLAGS = -lm  -lagbnp -lnblist  -L. 
OMPF=
SAVEF=

FFTWOMPLIB=-lfftw3_threads


# ===========  platform specific variables =====================
#  e.g. default compiler names, MPI wrapper name, etc

FC=mpif90
CC=gcc
LD=mpif90
FFLAGS :=  $(FFLAGS)   -ffree-line-length-0  -fallow-argument-mismatch
CFLAGS :=  $(CFLAGS)  -DFC_UNDERSCORES
F90FLAGS = -O3
LDFLAGS := $(LDFLAGS)   -lfftw3   -lfftw3_threads
OMPF:=$(OMPF) -fopenmp
SAVEF:=$(SAVEF) -fno-automatic

# ============================================================================
# NAMING conventions

# x.f           FORTRAN 
# x.f90         FORTRAN 90 code in free form
# [y.CPP.f]     cpp-preprocessed FORTRAN (created from y.f; y.CPP.f may be deleted)
# x.h           header (include) file
# y.H           cpp-pre-processable header file (will create y.h; y.H must not be
#			deleted)
# z.inc         cpp-pre-processable header file (will create z.h; z.inc must not be
#			deleted)


# ==== lists of objects, includes, etc. 

# objects of FORTRAN codes
OBJFOR = \
  add1h.o \
  add2h.o \
  add3h.o \
  add_char1.o \
  add_concta.o \
  add_energies.o \
  add_hill.o \
  addh2o.o \
  addh.o \
  add_int1.o \
  add_int2.o \
  addlnk.o \
  add_mapnl.o \
  add_real1.o \
  add_real2.o \
  add_solvent_coord.o \
  add_solvent_tpg.o \
  add_str.o \
  addtor.o \
  addtpg.o \
  anneal.o \
  appatm.o \
  appbkn.o \
  appbou.o \
  appct.o \
  appgrp.o \
  asng_xrmsw.o \
  assign_agbnp_types.o \
  assign_sasa_types.o \
  asstpe.o \
  atoms.o \
  atthd.o \
  bcnstp.o \
  bend.o \
  bldbox.o \
  bond.o \
  brent.o \
  brent_total.o \
  brot.o \
  calc_avg2_str.o \
  calc_avg_str.o \
  calc_avg_xrms.o \
  calc_gofr.o \
  calc_xrms.o \
  caldr.o \
  cc.o \
  change_coord_inner.o \
  change_frame.o \
  change_origin.o \
  change_step_sim.o \
  change_tors.o \
  chder_dyn.o \
  chebev_dyn.o \
  chebft_dyn.o \
  check_length_sim.o \
  check_read_columns.o \
  check_topology.o \
  check_zero_coord.o \
  chk14.o \
  chkewl.o \
  chlab.o \
  chmass.o \
  chnpr.o \
  close_files.o \
  clpcm.o \
  clt14.o \
  collision.o \
  comd0.o \
  comp_abmd_fdiss.o \
  comp_cell_center.o \
  comp_concth.o \
  comp_dip2.o \
  comp_dip.o \
  comp_dynamic_mat.o \
  comp_fcm.o \
  comp_forcep.o \
  comp_gmgp.o \
  comp_molmass.o \
  comp_nmol.o \
  comp_rmsq.o \
  comp_stress_bnd.o \
  comp_stress_conf.o \
  comp_stress_kinetic.o \
  comp_thermos_energy.o \
  comp_thermos_forces.o \
  comp_total_energy.o \
  comp_vcm.o \
  comp_vel_labframe.o \
  comp_work.o \
  coordinate_spline.o \
  coordinate_spline_init.o \
  copy_protl.o \
  correc.o \
  correc3x3.o \
  correc_etap.o \
  correc_exp_scale.o \
  correc_gmgp.o \
  correc_matr.o \
  correc_scale_2.o \
  correc_scale_co.o \
  correc_scale.o \
  correc_stress.o \
  covbd.o \
  covchg.o \
  covfod.o \
  covnb.o \
  covpr.o \
  cov_thermos.o \
  covuen.o \
  cpropar.o \
  crdrtf.o \
  cself.o \
  cself_alchemy.o \
  csort.o \
  dblas.o \
  dcd.o \
  dist_ijk.o \
  dlineq.o \
  do_dynamic_input.o \
  drive_analysis.o \
  drive_com_restraint.o \
  dumptp.o \
  duni.o \
  dupl_char1.o \
  dupl_concta.o \
  dupl_int1.o \
  dupl_int2.o \
  dupl_mapnl.o \
  dupl_real1.o \
  dupl_reali2.o \
  eigrs.o \
  erf_corr_cutoff.o \
  erfc_spline.o \
  erfc_spline_init.o \
  f1dim_adjust_bonds.o \
  f1dim_der.o \
  fft3-public.o \
  fftpack.o \
  find_diss_mol.o \
  find_length_run.o \
  find_thermos.o \
  fndfmt.o \
  fndgrp.o \
  fpbond_adjust_bonds.o \
  gen_abmd_kvect.o \
  generate_rnd_rot.o \
  generate_slv.o \
  getchr.o \
  get_prot_cnstr.o \
  get_sgr.o \
  get_spectra_vacf.o \
  get_type_slv.o \
  get_velocities.o \
  hapacp.o \
  hapdon.o \
  igmap.o \
  indexx.o \
  initialize_lambda.o \
  initialize_work.o \
  inicmp.o \
  insert_solute.o \
  int_corr_erf_spline.o \
  interpol_dyna.o \
  jartemp.o \
  lc_index.o \
  link.o \
  linmin_adjust_bonds.o \
  linmin_total.o \
  low_up.o \
  lsqqtn.o \
  mappa.o \
  matinv.o \
  meta.o\
  meta_forces.o \
  meta_pot.o \
  meta_rc.o \
  min_pack.o \
  minv.o \
  mnbrakf.o \
  mnbrakf_total.o \
  mts_furier_alchemy.o \
  mts_plot_fragm.o \
  mts_plotp.o \
  mts_umbrl.o \
  near0.o \
  normal.o \
  openf.o \
  open_file_dump.o \
  parse.o \
  parse_numbers.o \
  pbend.o \
  pbond.o \
  pckben.o \
  pckfrg.o \
  pckitor.o \
  pcktor.o \
  perc_bondh.o \
  phapp.o \
  pitors.o \
  plot_center.o \
  plotd.o \
  pntres.o \
  pptors.o \
  print_alchemic_work.o \
  print_cm.o \
  print_title_analysis.o \
  prtat.o \
  prtba.o \
  prtbnd.o \
  prtcn.o \
  prtfrc.o \
  prtite.o \
  prtite_min.o \
  prtit.o \
  prtpt.o \
  prtsq.o \
  ranf.o \
  ranvel.o \
  rattle_correc_co.o \
  rattle_verlet_co.o \
  rdaprm.o \
  rdatpg.o \
  rdcmac.o \
  rdrgd.o \
  read_analysis.o \
  read_confc_columns.o \
  read_confc_rows.o \
  read_integrator.o \
  read_meta.o \
  read_parameters.o \
  read_potential.o \
  read_properties.o \
  read_run.o \
  read_simulation.o \
  read_solvent.o \
  read_space_group.o \
  readtp.o \
  recpr.o \
  redefine_slv.o \
  redmss.o \
  remove_mv.o \
  repl_coord.o \
  repl_tpg.o \
  rigfit.o \
  rmbend.o \
  rmpbond.o \
  rotb.o \
  sasamod.o \
  sbend.o \
  sbond.o \
  scale_charges.o \
  scale_masses.o \
  sctab.o \
  search_clsth.o \
  serchd.o \
  set_const_co.o \
  set_constxy_co.o \
  set_hbonds_masks.o \
  set_ss_array.o \
  set_tempp.o \
  set_tempt.o \
  setuns.o \
  setup_box.o \
  setup_skin_shell.o \
  spec14.o \
  starting_verlet.o \
  steer_along_path.o \
  alchemic.o\
  start.o \
  steer_dynamics.o \
  steer_print.o \
  steer_variables.o \
  strblk.o \
  time_correlation.o \
  tors.o \
  tr_inbox.o \
  twofft.o \
  unit.o \
  update.o \
  up_low.o \
  uscrpl.o \
  usort.o \
  verify_variables.o \
  verlet.o \
  verlet_free_eta.o \
  verlet_free.o \
  windows.o \
  write_bends.o \
  write_bonds.o \
  write_confc.o \
  write_diffusion.o \
  write_fragm_dist.o \
  write_gofrp.o \
  write_gofrw.o \
  write_gyr.o \
  write_hbonds.o \
  write_pot_bond.o \
  write_pot_nbond.o \
  write_rms.o \
  write_tors.o \
  write_vacf_phi.o \
  write_xrms_atm.o \
  write_xrms.o \
  xfit.o \
  zero0.o \
  zero3x3.o \
  zeroa.o \
  zero.o \
  zero_gofr.o \
# objects that must be precompiled with cpp
OBJPRE= \
  fourier.CPP.o\
  furippa.CPP.o \
  furipp_alchemy.CPP.o \
  giofar.CPP.o \
  parst.CPP.o \
  spme.CPP.o \
  pme_fftw3.CPP.o \
  omp_integr.o
# objects specific of MPI version
OBJPARA= \
  read_rem.o \
  orac_mpi.CPP.o \
# objects that need to be precompiled (SERIAL/MPI/FFTW/OMP flavours) 
OBJPARAPRE= \
  adjust_bonds.CPP.o \
  alchemy_charge_grid.CPP.o \
  blkdta.CPP.o         \
  do_pmesh_alchemy.CPP.o \
  dumprs.CPP.o         \
  ferrf.CPP.o \
  ferrf_alchemy.CPP.o \
  ferrf_alchemy_add.CPP.o \
  f1dim_total.CPP.o \
  fft_pme.CPP.o \
  fft_pme_alchemy.CPP.o \
  fftw3_pme_init.CPP.o \
  gauss.CPP.o \
  get_total_energy.CPP.o \
  force_alchemy.CPP.o \
  force_alchemy_shifted.CPP.o \
  force.CPP.o \
  force_REM.CPP.o \
  fpbond.CPP.o \
  fpbend.CPP.o \
  fptors.CPP.o \
  fptors_rem.CPP.o \
  fitors.CPP.o \
  fitors_rem.CPP.o \
  fnb14.CPP.o \
  fnb14_alchemy.CPP.o \
  fnb14_rem.CPP.o \
  fnbgrp.CPP.o \
  fnbgrp_alchemy.CPP.o \
  join.CPP.o \
  kinetic.CPP.o \
  lc_list.CPP.o \
  mts_forces.CPP.o \
  mts_intra_n0.CPP.o \
  mts_intra_n1.CPP.o \
  mts_furier.CPP.o \
  mts_test.CPP.o       \
  mtsmd.CPP.o          \
  orac.CPP.o           \
  plotc.CPP.o          \
  prtacc.CPP.o          \
  rattle_correc.CPP.o \
  rattle_verlet.CPP.o \
  read_inout.CPP.o     \
  read_input.CPP.o     \
  read_setup.CPP.o \
  read_solute.CPP.o \
  readco.CPP.o         \
  readrs.CPP.o         \
  rem.CPP.o            \
  run_minimize.CPP.o \
  sasa.CPP.o \
  scalar_sum_alchemy.CPP.o \
  scalar_sum.CPP.o \
  sge.CPP.o            \
  timer.CPP.o \
  verify_input.CPP.o   \
  vi_potential.CPP.o \
  wrenc.CPP.o \
  xerror.CPP.o

#-------------------------------------------------------
#define OBJ with !$OMP directives and NOT -save option
#------------------------------------------------------

OBJF_OMP= \
  appbou.o \
  change_frame.o \
  comp_stress_kinetic.o \
  comp_thermos_force.o \
  correc.o \
  correc_gmgp.o \
  correc_scale.o \
  correc_scale_2.o \
  verlet.o \
  ferrf.CPP.o \
  ferrf_alchemy.CPP.o \
  ferrf_alchemy_add.CPP.o \
  fitors.CPP.o \
  fitors_rem.CPP.o \
  fnb14.CPP.o \
  fnb14_alchemy.CPP.o \
  fnb14_rem.CPP.o \
  fnbgrp.CPP.o \
  fnbgrp_alchemy.CPP.o \
  force.CPP.o \
  force_alchemy.CPP.o \
  force_alchemy_shifted.CPP.o \
  force_REM.CPP.o \
  fpbond.CPP.o \
  fpbend.CPP.o \
  fptors.CPP.o \
  fptors_rem.CPP.o \
  kinetic.CPP.o \
  lc_list.CPP.o \
  mts_intra_n0.CPP.o \
  mts_intra_n1.CPP.o \
  rattle_correc.CPP.o \
  rattle_verlet.CPP.o \
  read_input.CPP.o \
  scalar_sum.CPP.o \
  timer.CPP.o \
  verlet_free.o 

OBJF90_OMP= \
  alchemy_charge_grid.CPP.o \
  do_pmesh_alchemy.CPP.o \
  pme_fftw3.CPP.o \
  omp_integr.o \
  sasa.CPP.o \
  scalar_sum_alchemy.CPP.o \

#-------------------------------------------------------
#define OBJ with !$OMP directives AND -save option
#------------------------------------------------------

OBJF_OMP_SAVE= \
  fft_pme.CPP.o \
  fft_pme_alchemy.CPP.o \
  mtsmd.CPP.o\
  run_minimize.CPP.o\
  orac.CPP.o

OBJF90_OMP_SAVE= \

#-------------------------------------------------------
#define OBJ with -save option only
#------------------------------------------------------

OBJF_SAVE= \
  mts_furier.CPP.o \
  jartemp.o \

OBJF90_SAVE= \

#-------------------------------------------------------
#       sets rules for compiler FLAGS accordingly.
#------------------------------------------------------

# compile with OMP only 
  $(OBJF_OMP) : %.o : %.f;  $(FC) -c $< $(FFLAGS) $(OMPF) -o $@
  $(OBJF90_OMP) : %.o : %.f90;  $(FC) -c $< $(FFLAGS) $(F90FLAGS) $(OMPF) -o $@
# compile with -save only 
  $(OBJF_SAVE) : %.o : %.f;  $(FC) -c $< $(FFLAGS) $(SAVEF) -o $@
  $(OBJF90_SAVE) : %.o : %.f90;  $(FC) -c $< $(FFLAGS) $(SAVEF) -o $@
# compile with -save AND OMP 
  $(OBJF_OMP_SAVE) : %.o : %.f;  $(FC) -c $< $(FFLAGS) $(SAVEF) $(OMPF) -o $@
  $(OBJF90_OMP_SAVE) : %.o : %.f90;  $(FC) -c $< $(FFLAGS) $(SAVEF) $(OMPF) -o $@


# objects for the version including agbnp
OBJ_AGBNP = agbnp.o

# libraries for the version including agbnp
LIB_AGBNP_1 = libagbnp.a
LIB_AGBNP_2 = libnblist.a

# precious header files
HEADER= \
  config.H \
  pressure.H

# .h/.f always included
INCLUDEARCH= \
  orac.inc \
  pbc.inc \
  parameters.inc lc_list.inc \
  iobuffer.inc \
  mtsmd_avg_inc.f \
  mtsmd_dump_inc.f
# list of source codes
LISTFORC = \
  $(patsubst %.o, %.f, $(OBJFOR)) \
  $(patsubst %.o, %.c, $(OBJC)) \
  $(patsubst %.o, %.f, $(OBJPARA)) \
  $(patsubst %.o, %.f, $(OBJOMP)) \
  $(patsubst %.CPP.o, %.f, $(OBJPRE)) \
  $(patsubst %.CPP.o, %.f, $(OBJPARAPRE)) 

# complete list of source files
LIST_NO_MAKEFILE = \
  configure \
  config.h_dist \
  $(INCLUDEARCH) \
  $(INCLUDEPARE) \
  $(HEADER) \
  $(LISTFORC)

LIST = \
  Makefile \
  $(LIST_NO_MAKEFILE)


ARFLAGS = rs
AGBNP_OBJ = libagbnp.o libagbnp_pbc.o


# ===========================================================================
#  MACROS
# ===========================================================================


libnblist.a:  libnblist.o
	$(AR) $(ARFLAGS) libnblist.a libnblist.o

libnblist.o: libnblist.c libnblist.h
	$(CC) -c $(CFLAGS) $<


clean_nbl: 
	rm libnblist.a libnblist.o



libagbnp.a:  $(AGBNP_OBJ) libnblist.a
	$(AR) $(ARFLAGS) libagbnp.a $(AGBNP_OBJ)

libagbnp.o: libagbnp.c agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -c -DUSE_MMNBLIST $<

libagbnp_pbc.o: libagbnp.c agbnp.h agbnp_private.h
	$(CC) $(CFLAGS) -c -DUSE_PBC -DUSE_ALT -DUSE_MMNBLIST -o libagbnp_pbc.o $<

clean_agbnp: 

	rm libagbnp.*

clean_lib: clean_agbnp clean_nbl


# ===========================================================================
# === USER-LEVEL TARGETS ====================================================
# ===========================================================================

# since target is a directory name, define it as "phony" to `make' it anyway
# default: compile and link in subdir "default", using defaults for FC etc

.PHONY: GNU-FFTW-OMP-MPI
GNU-FFTW-OMP-MPI:  
	@$(RM) *.CPP.o *.a *.mod *.CPP.f*  # remove temp file in main DIR
	-mkdir GNU-FFTW-OMP-MPI
	$(MAKE) -C GNU-FFTW-OMP-MPI MAKEFILES="../Makefile" \
	"CPPARG = $(CPPARG)   -DFFTW -D_OMP_ -D_MPI_" \
        "FC = mpif90" \
        "FFLAGS :=  $(FFLAGS)"\
        "CFLAGS :=  $(CFLAGS)"\
         orac


# ===============================================================================
# === MAIN COMPILATION RULE                    ==================================
# ===============================================================================

# link the big exe program (with many source files)

ifeq ($(word 1,$(FC)), mpif90)
  MPILIST=$(OBJPARA) 
else
  ifeq ($(word 1,$(FC)), mpiifort)
     MPILIST=$(OBJPARA) 
  else
     MPILIST=
  endif
endif



#       link the big exe program (for MPI COMPILATION
orac:  $(OBJFOR) $(OBJPRE) $(OBJPARAPRE) $(OBJ_AGBNP) $(LIB_AGBNP_1) $(LIB_AGBNP_2) $(MPILIST)
	$(FC) -o orac $(OBJFOR) $(OBJPRE) $(OBJPARAPRE) $(OBJ_AGBNP) $(MPILIST)  $(LDFLAGS)  $(OMPF) 

# === CANCEL RCS IMPLICIT RULE                 ==================================
# ===============================================================================

# this disables searching of prerequisites in RCS/ subdir by default
# to enable RCS search, type `make RCS=1 [target]'
ifndef RCS
%: RCS/%,v
endif

# ===============================================================================
# === ADD .f90 IMPLICIT RULE                   ==================================
# ===============================================================================

#
# N.B. Always use .f90 suffix for a FORTRAN90 code in FREE FORM
#

# don't use IMPLICIT rules like the following, as Make is confused by
#   the existence of .mod files, and will apply the Modula-2 implicit
#   rule when you don't want it to

#.SUFFIXES: .f90
#.f90.o:
#	$(FC) $(FFLAGS) $(F90FLAGS) -c $<


# note that gfortran won't touch "%.mod" after recompiling 
#      if it's not changed; 
#      to force unconditional touch, uncomment last line
# see http://lagrange.mechse.illinois.edu/mwest/f90_mod_deps/
%.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -c $<
#	touch $(patsubst %.f90,%.mod,$<)

# ===============================================================================
# === OBJECT DEPENDENCIES                      ==================================
# ===============================================================================

# use a special compile flag for this file
#recpr.o: recpr.f90
#ifdef FFLAGS1
#	$(FC) $(FFLAGS1) -c -o $@ $<
#else	
#	$(FC) $(FFLAGS) -c -o $@ $<
#endif

# .CPP.o: .CPP.f  is implicit 
# .CPP.o: .CPP.f90  is NOT implicit 
add_hill.o: add_hill.f unit.o meta.o
add_solvent_tpg.o: add_solvent_tpg.f cpropar.o parst.CPP.o
adjust_bonds.CPP.o: adjust_bonds.f parst.CPP.o
appbou.o: appbou.f omp_integr.o
assign_sasa_types.o: assign_sasa_types.f unit.o
bcnstp.o: bcnstp.f pbc.h
bldbox.o: bldbox.f cpropar.o parst.CPP.o
ifeq ($(findstring MPI,$(CPPARG)), MPI)
blkdta.CPP.o: blkdta.f90 cpropar.o parst.CPP.o unit.o rem.CPP.o sge.CPP.o meta.o dcd.o orac_mpi.CPP.o
else
blkdta.CPP.o: blkdta.f90 cpropar.o parst.CPP.o unit.o rem.CPP.o sge.CPP.o meta.o dcd.o
endif
brot.o: brot.f unit.o
calc_avg2_str.o: calc_avg2_str.f parst.CPP.o
calc_avg_str.o: calc_avg_str.f parst.CPP.o
calc_gofr.o: calc_gofr.f parst.CPP.o pbc.h
calc_xrms.o: calc_xrms.f parst.CPP.o
change_tors.o: change_tors.f parameters.h parst.CPP.o
chkewl.o: chkewl.f fourier.CPP.o unit.o
close_files.o: close_files.f unit.o
collision.o: collision.f parst.CPP.o unit.o
comp_dip.o: comp_dip.f pbc.h
comp_dynamic_mat.o: comp_dynamic_mat.f cpropar.o parst.CPP.o unit.o
comp_thermos_forces.o: comp_thermos_forces.f unit.o
comp_total_energy.o: comp_total_energy.f cpropar.o parst.CPP.o unit.o
correc.o: correc.f omp_integr.o
correc_stress.o: correc_stress.f unit.o
covbd.o: covbd.f unit.o
covchg.o: covchg.f unit.o
covfod.o: covfod.f unit.o
covnb.o: covnb.f unit.o
covpr.o: covpr.f unit.o
cov_thermos.o: cov_thermos.f unit.o
covuen.o: covuen.f unit.o
cpropar.o: cpropar.f90 parst.CPP.o
crdrtf.o: crdrtf.f cpropar.o parst.CPP.o unit.o
cself.o: cself.f unit.o
cself_alchemy.o: cself_alchemy.f unit.o
dcd.o: dcd.f90
do_dynamic_input.o: do_dynamic_input.f cpropar.o parst.CPP.o
do_pmesh_alchemy.CPP.o:   do_pmesh_alchemy.f90 spme.CPP.o 
drive_analysis.o: drive_analysis.f cpropar.o fourier.CPP.o giofar.CPP.o iobuffer.h lc_list.h parst.CPP.o spme.CPP.o unit.o
dumprs.CPP.o: dumprs.f rem.CPP.o sge.CPP.o
dumptp.o: dumptp.f cpropar.o parst.CPP.o
erf_corr_cutoff.o: erf_corr_cutoff.f unit.o
f1dim_adjust_bonds.o: f1dim_adjust_bonds.f cpropar.o parst.CPP.o
f1dim_der.o: f1dim_der.f cpropar.o parst.CPP.o
f1dim_total.CPP.o: f1dim_total.f cpropar.o parst.CPP.o
ferrf.CPP.o: ferrf.f unit.o omp_integr.o
ferrf_alchemy.CPP.o: ferrf_alchemy.f omp_integr.o unit.o
ferrf_alchemy_add.CPP.o: ferrf_alchemy_add.f unit.o omp_integr.o
fft_pme.CPP.o: fft_pme.f spme.CPP.o
fftw3_pme_init.o:  fftw3_pme_init.f90  omp_integr.o
# fft_pme_init.CPP.o:  fft_pme_init.CPP.f spme.CPP.o
fft_pme_alchemy.CPP.o:  fft_pme_alchemy.f spme.CPP.o
fft_pme_alchemy.o: fft_pme_alchemy.f spme.CPP.o
find_length_run.o: find_length_run.f unit.o
fnb14.CPP.o: fnb14.f unit.o  omp_integr.o
fnb14_alchemy.CPP.o: fnb14_alchemy.f unit.o  omp_integr.o
fnb14_rem.CPP.o: fnb14_rem.f unit.o rem.CPP.o  omp_integr.o
force_alchemy.CPP.o: force_alchemy.f unit.o parst.CPP.o  omp_integr.o
force_alchemy_shifted.CPP.o: force_alchemy_shifted.f unit.o parst.CPP.o  omp_integr.o
force.CPP.o: force.f unit.o parst.CPP.o  cpropar.o
force_REM.CPP.o: force_REM.f unit.o parst.CPP.o omp_integr.o
fnbgrp.CPP.o: fnbgrp.f pbc.h unit.o omp_integr.o
fnbgrp_alchemy.CPP.o: fnbgrp_alchemy.f pbc.h unit.o omp_integr.o
fpbond.CPP.o: fpbond.f omp_integr.o
fpbend.CPP.o: fpbend.f omp_integr.o
fptors.CPP.o: fptors.f omp_integr.o
fptors_rem.CPP.o: fptors_rem.f rem.f90 rem.CPP.o  omp_integr.o
fitors.CPP.o: fitors.f unit.o omp_integr.o
fitors_rem.CPP.o: fitors_rem.f rem.f90 rem.CPP.o  omp_integr.o
ifeq ($(findstring MPI,$(CPPARG)), MPI)
  gauss.CPP.o: gauss.f90  rem.CPP.o unit.o orac_mpi.CPP.o
else
  gauss.CPP.o: gauss.f90  unit.o rem.CPP.o 
endif
etchr.o: getchr.f unit.o
get_total_energy.CPP.o: get_total_energy.f cpropar.o parst.CPP.o unit.o sasamod.o
initialize_work.o: initialize_work.f cpropar.o parst.CPP.o spme.CPP.o
insert_solute.o: insert_solute.f pbc.h
join.CPP.o: join.f parameters.h parst.CPP.o unit.o
kinetic.CPP.o: kinetic.f parst.CPP.o unit.o omp_integr.o
lc_list.CPP.o: lc_list.f lc_list.h parst.CPP.o pbc.h omp_integr.o
meta_forces.CPP.o: meta_forces.f unit.o meta.o
meta_pot.o: meta_pot.f unit.o meta.o
meta_rc.o: meta_rc.f pbc.h
mts_forces.o: mts_forces.f cpropar.o parst.CPP.o unit.o rem.CPP.o omp_integr.o
mts_furier.CPP.o: mts_furier.f cpropar.o parst.CPP.o spme.CPP.o omp_integr.o
mts_furier_alchemy.o: mts_furier_alchemy.f cpropar.o parst.CPP.o spme.CPP.o omp_integr.o
ifeq ($(findstring MPI,$(CPPARG)), MPI)
  mts_intra_n0.CPP.o: mts_intra_n0.f cpropar.o parst.CPP.o rem.CPP.o pbc.h  orac_mpi.CPP.o unit.o omp_integr.o
  mts_intra_n1.CPP.o: mts_intra_n1.f cpropar.o parst.CPP.o rem.CPP.o meta.o orac_mpi.CPP.o unit.o omp_integr.o
else
  mts_intra_n0.CPP.o: mts_intra_n0.f cpropar.o parst.CPP.o rem.CPP.o pbc.h  omp_integr.o
  mts_intra_n1.CPP.o: mts_intra_n1.f cpropar.o parst.CPP.o rem.CPP.o meta.o omp_integr.o
endif
mts_test.CPP.o: mts_test.f sge.CPP.o
mtsmd.CPP.o: mtsmd.f cpropar.o fourier.CPP.o giofar.CPP.o lc_list.h mtsmd_avg_inc.CPP.f mtsmd_dump_inc.CPP.f parst.CPP.o spme.CPP.o unit.o rem.CPP.o sge.CPP.o meta.o dcd.o omp_integr.o
mts_plot_fragm.o: mts_plot_fragm.f unit.o
mts_plotp.o: mts_plotp.f unit.o
mts_umbrl.o: mts_umbrl.f parst.CPP.o pbc.h
omp_integr.o: omp_integr.f90 parst.CPP.o 
open_file_dump.o: open_file_dump.f unit.o
ifeq ($(findstring MPI,$(CPPARG)), MPI)
  orac.CPP.o: orac.f cpropar.o orac.h parst.CPP.o unit.o rem.CPP.o  orac_mpi.CPP.o omp_integr.o
else 
  orac.CPP.o: orac.f cpropar.o orac.h parst.CPP.o unit.o rem.CPP.o  omp_integr.o
endif
orac_mpi.CPP.o: orac_mpi.f90 rem.CPP.o cpropar.o omp_integr.o
pckben.o: pckben.f parst.CPP.o
pckitor.o: pckitor.f parst.CPP.o
pcktor.o: pcktor.f parst.CPP.o
pitors.o: pitors.f parst.CPP.o
plot_center.o: plot_center.f unit.o
plotc.CPP.o: plotc.f unit.o rem.CPP.o sge.CPP.o
pme_fftw3.CPP.o:  pme_fftw3.f90
pptors.o: pptors.f parst.CPP.o
print_title_analysis.o: print_title_analysis.f cpropar.o parst.CPP.o unit.o
prtacc.CPP.o: prtacc.f cpropar.o parst.CPP.o unit.o
prtat.o: prtat.f unit.o
prtba.o: prtba.f pbc.h unit.o
prtbnd.o: prtbnd.f pbc.h unit.o
prtcn.o: prtcn.f pbc.h unit.o
prtite_min.o: prtite_min.f cpropar.o parst.CPP.o unit.o
prtite.o: prtite.f cpropar.o parst.CPP.o unit.o
prtit.o: prtit.f pbc.h unit.o
prtpt.o: prtpt.f pbc.h unit.o
prtsq.o: prtsq.f parst.CPP.o unit.o
ranvel.o: ranvel.f parst.CPP.o unit.o
rattle_correc.CPP.o: rattle_correc.f parst.CPP.o omp_integr.o
rattle_verlet.CPP.o: rattle_verlet.f parst.CPP.o omp_integr.o
rdcmac.o: rdcmac.f parst.CPP.o cpropar.o rem.CPP.o
read_analysis.o: read_analysis.f cpropar.o parst.CPP.o unit.o
read_confc_columns.o: read_confc_columns.f unit.o
read_confc_rows.o: read_confc_rows.f unit.o
read_inout.CPP.o: read_inout.f cpropar.o parst.CPP.o unit.o rem.CPP.o dcd.o omp_integr.o
read_input.CPP.o: read_input.f cpropar.o parst.CPP.o unit.o rem.CPP.o sge.CPP.o omp_integr.o
read_integrator.o: read_integrator.f cpropar.o parst.CPP.o unit.o
read_meta.o: read_meta.f parst.CPP.o cpropar.o unit.o meta.o
read_parameters.o: read_parameters.f parameters.h parst.CPP.o
read_potential.o: read_potential.f cpropar.o parst.CPP.o unit.o
read_properties.o: read_properties.f cpropar.o parst.CPP.o unit.o
read_rem.o: read_rem.f cpropar.o parst.CPP.o unit.o rem.CPP.o
read_run.o: read_run.f cpropar.o parst.CPP.o unit.o
read_simulation.o: read_simulation.f cpropar.o parst.CPP.o unit.o
ifeq ($(findstring MPI,$(CPPARG)), MPI)
  read_setup.CPP.o: read_setup.f cpropar.o parst.CPP.o unit.o orac_mpi.CPP.o
  read_solute.CPP.o: read_solute.f cpropar.o parst.CPP.o unit.o orac_mpi.CPP.o
else
  read_setup.CPP.o: read_setup.f cpropar.o parst.CPP.o unit.o 
  read_solute.CPP.o: read_solute.f cpropar.o parst.CPP.o unit.o
endif
read_solvent.o: read_solvent.f cpropar.o parst.CPP.o unit.o
readco.CPP.o: readco.f unit.o rem.CPP.o sge.CPP.o
readrs.CPP.o: readrs.f rem.CPP.o sge.CPP.o
readtp.o: readtp.f cpropar.o parst.CPP.o
rem.CPP.o: rem.f90 unit.o dcd.o
repl_coord.o: repl_coord.f unit.o
repl_tpg.o: repl_tpg.f cpropar.o parst.CPP.o
rmbend.o: rmbend.f parst.CPP.o
rotb.o: rotb.f unit.o
run_minimize.CPP.o: run_minimize.f cpropar.o lc_list.h parst.CPP.o spme.CPP.o unit.o sasamod.o
sasamod.o: sasamod.f90 parst.CPP.o
set_tempp.o: set_tempp.f unit.o
set_tempt.o: set_tempt.f unit.o
setuns.o: setuns.f unit.o
setup_box.o: setup_box.f cpropar.o parst.CPP.o unit.o
sge.CPP.o: sge.f90 cpropar.o unit.o rem.CPP.o
start.o: start.f cpropar.o parst.CPP.o unit.o
tr_inbox.o: tr_inbox.f pbc.h
unit.o: unit.f90 
update.o: update.f parst.CPP.o pbc.h unit.o
verify_input.CPP.o: verify_input.f cpropar.o parst.CPP.o spme.CPP.o unit.o rem.CPP.o sge.CPP.o dcd.o
verify_variables.o: verify_variables.f cpropar.o parst.CPP.o unit.o
verlet.o: verlet.f omp_integr.o
vi_potential.CPP.o: vi_potential.f cpropar.o parst.CPP.o pbc.h spme.CPP.o unit.o
write_bends.o: write_bends.f parst.CPP.o
write_bonds.o: write_bonds.f parst.CPP.o
write_confc.o: write_confc.f unit.o
write_diffusion.o: write_diffusion.f unit.o
write_fragm_dist.o: write_fragm_dist.f pbc.h
write_gofrp.o: write_gofrp.f unit.o
write_gofrw.o: write_gofrw.f unit.o
write_gyr.o: write_gyr.f parst.CPP.o
write_hbonds.o: write_hbonds.f parst.CPP.o pbc.h
write_pot_bond.o: write_pot_bond.f parst.CPP.o unit.o
write_pot_nbond.o: write_pot_nbond.f parst.CPP.o unit.o
write_tors.o: write_tors.f parst.CPP.o
write_vacf_phi.o: write_vacf_phi.f unit.o
write_xrms.o: write_xrms.f parst.CPP.o


agbnp.o: agbnp.c $(LIB_AGBNP_1) $(LIB_AGBNP_2)

# Preprocessable .f dependencies and rules 

adjust_bonds.CPP.f: adjust_bonds.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
alchemy_charge_grid.CPP.f90: alchemy_charge_grid.f90 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
blkdta.CPP.f90: blkdta.f90 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
do_pmesh_alchemy.CPP.f90: do_pmesh_alchemy.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
dumprs.CPP.f: dumprs.f	
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
f1dim_total.CPP.f: f1dim_total.f	
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fft_pme.CPP.f: fft_pme.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fft_pme_alchemy.CPP.f: fft_pme_alchemy.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fftw3_pme_init.CPP.f90: fftw3_pme_init.f90 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
force_alchemy.CPP.f: force_alchemy.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
force_alchemy_shifted.CPP.f: force_alchemy_shifted.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
force.CPP.f: force.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
force_REM.CPP.f: force_REM.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fourier.CPP.f90: fourier.f90 config.H
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fnb14.CPP.f: fnb14.f unit.o 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fnb14_alchemy.CPP.f: fnb14_alchemy.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fnb14_rem.CPP.f: fnb14_rem.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fnbgrp.CPP.f: fnbgrp.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fnbgrp_alchemy.CPP.f: fnbgrp_alchemy.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fpbond.CPP.f: fpbond.f  
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fpbend.CPP.f: fpbend.f  
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
ferrf.CPP.f: ferrf.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
ferrf_alchemy.CPP.f: ferrf_alchemy.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
ferrf_alchemy_add.CPP.f: ferrf_alchemy_add.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fitors.CPP.f: fitors.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fitors_rem.CPP.f: fitors_rem.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fptors.CPP.f: fptors.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
fptors_rem.CPP.f: fptors_rem.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
furippa.CPP.f: furippa.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
furipp_alchemy.CPP.f: furipp_alchemy.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
get_total_energy.CPP.f: get_total_energy.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
giofar.CPP.f90: giofar.f90 config.H
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
gauss.CPP.f90: gauss.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
join.CPP.f: join.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
kinetic.CPP.f: kinetic.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
lc_list.CPP.f: lc_list.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mts_forces.CPP.f: mts_forces.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mts_intra_n0.CPP.f: mts_intra_n0.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mts_intra_n1.CPP.f: mts_intra_n1.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mts_furier.CPP.f: mts_furier.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mts_test.CPP.f: mts_test.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mtsmd.CPP.f: mtsmd.f mtsmd_avg_inc.f mtsmd_dump_inc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mtsmd_avg_inc.CPP.f: mtsmd_avg_inc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
mtsmd_dump_inc.CPP.f: mtsmd_dump_inc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
orac.CPP.f: orac.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
orac_mpi.CPP.f90: orac_mpi.f90 rem.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
parst.CPP.f90: parst.f90 config.H
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
plotc.CPP.f: plotc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
#pme_code.CPP.f: pme_code.f
#	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
pme_fftw3.CPP.f90: pme_fftw3.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
prtacc.CPP.f: prtacc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
read_inout.CPP.f: read_inout.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
rattle_correc.CPP.f: rattle_correc.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
rattle_verlet.CPP.f: rattle_verlet.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
read_input.CPP.f: read_input.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
read_setup.CPP.f: read_setup.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
read_solute.CPP.f: read_solute.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
readco.CPP.f: readco.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
readrs.CPP.f: readrs.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
rem.CPP.f90: rem.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
run_minimize.CPP.f: run_minimize.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
sasa.CPP.f90: sasa.f90 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
scalar_sum.CPP.f: scalar_sum.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
scalar_sum_alchemy.CPP.f90: scalar_sum_alchemy.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
sge.CPP.f90: sge.f90
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
spme.CPP.f90: spme.f90 config.H
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
timer.CPP.f: timer.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
verify_input.CPP.f: verify_input.f pressure.H
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
vi_potential.CPP.f: vi_potential.f
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
wrenc.CPP.f: wrenc.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@
xerror.CPP.f: xerror.f 
	$(CPP) $(CPPARG) $< | awk '!/^ *$$/ && !/^#/' > $@


# Preprocessable includes dependencies and rules

iobuffer.h: iobuffer.inc config.H
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@
lc_list.h: lc_list.inc config.H
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@
orac.h: orac.inc
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@
parameters.h: parameters.inc
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@
pbc.h: pbc.inc config.H
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@
voronoi.CPP.f90: voronoi.f90 config.H
	$(CPP) $< | awk '!/^ *$$/ && !/^#/' > $@

# ===========================================================================
# === DEVELOPMENT AND SERVICE ===============================================
# ===========================================================================


# clean: remove useless files, but keep the big executable
clean:
	@$(RM) core ?*[~#] *.o *.a *.mod *.CPP.f* *.h temp

# "clean" each directory; to clean only one, type: `make DIRS=<dir> dirclean'
dirclean: 
	for dir in $(DIRS); do \
	   $(MAKE) MAKEFILES=../Makefile -C $$dir clean ;  \
	done

# remove also the big executables and directories
realclean: clean
	@$(RM) -rf orac $(DIRS)

distclean: realclean

testorac:
	@echo "FC =" $(FC)
	@echo "FFLAGS =" $(FFLAGS)
	@echo "CPPARG =" $(CPPARG)
	@echo "OMPF =" $(OMPF)
	@echo "OMPF =" $(SAVEF)

test:
	$(MAKE) "CPPARG = -D_MPI_ $(CPPARG)" testorac
