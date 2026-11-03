      SUBROUTINE mtsmd(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm,ypcm,zpcm
     &     ,ndim)
************************************************************************
*                                                                      *
*     MTSMD is the driver of the MD run when multiple time scales      *
*     methods are used. The integration algorithm is r-RESPA           *
*     of order o(dt^3) (Tuckermann et al. JCP 97 1990 (1992))          *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*     New constant pressure and temperature algorithms implemented     *
*     by Massimo Marchi                                                *
*                                                                      *
*______________________________________________________________________*
*                                                                      *
*  MTSMD externals:                                                    *
*     adjust_bonds anneal appbou asng_xrmsw                            *  
*     bcnstp calc_avg_str calc_avg_xrms calc_gofr calc_xrms            *
*     change_coord_inner change_frame change_origin check_topology     *
*     chkewl collision comp_abmd_fdiss comp_dip comp_fcm               *
*     comp_forcep comp_molmass comp_stress_conf                        *
*     comp_stress_kinetic comp_thermos_energy comp_thermos_forces      *
*     comp_vcm comp_vel_labframe copy_protl correc                     *
*     correc_etap correc_exp_scale correc_matr correc_stress cov_thermo*
*     cself daxpy dcopy dscal dumprs                                   *
*     erfc_spline_init erf_corr_cutoff ferrf fft_pme_init              *
*     find_diss_mol find_thermos fndgrp gen_abmd_kvect                 *
*     get_prot_cnstr get_type_slv icopy                                *
*     inicmp int_corr_erf_spline kinetic lc_index lc_list              *
*     matinv mts_forces mts_furier mts_intra_n0   mts_intra_n1         *
*     mts_plotp mts_plot_fragm mts_test plotc plotd                    *
*     plot_center prtacc prtat prtba prtbnd                            *
*     prtcn prtfrc prtit prtite prtpt prtsq                            *
*     ranvel rattle_correc rattle_correc_co rattle_verlet              *
*     rattle_verlet_co readrs setup_skin_shell set_const_co            *
*     set_tempp set_tempt starting_verlet timer trans_center           *
*     tr_inbox verlet verlet_free verlet_free_eta write_bends          *
*     write_bonds write_confc write_gofrp write_gofrw write_pot_bond   *
*     write_pot_nbond write_tors xerror zero zero0                     *
*     zero3x3 zeroa zero_gofr
* 								       *
************************************************************************


*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar
      use omp_integr
      use fourier
      use giofar
      use spme
      use meta
      use rem
      use sge, only: sge_drive, sge_print_titles, sge_init,
     &     sim_tempering, sim_tempering_sge,always_accept,ncycle
     &     ,ncycle_reject,use_saved,go_save,target_replica,sge_ts
     &     ,nfreqsge,print_always
      use dcd, only: ndcd,print_dcd
#ifdef _MPI_
      use orac_mpi
#endif
#ifdef _OMP_
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#endif
      IMPLICIT none
      
      INTEGER ma,mb,mh,ms,mf,memp,t1,mbs,numpr
      PARAMETER (ma=tsites,mb=tsitep,mh=elsiz,ms=secsiz,memp=1
     &     ,mf=mb*4,t1=types*(types+1)/2)

      PARAMETER (mbs=2*mb, numpr=npm)

*----------------------- ARGUMENTS -------------------------------------

#ifdef _OMP_
      INTEGER mapnl(m8t,nthr),ndim,cnstp_threads(nthr)
#else
      INTEGER mapnl(*),ndim
#endif
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),eta(*),xpcm(*)
     &     ,ypcm(*),zpcm(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      include 'lc_list.h'
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif

      
*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      EXTERNAL  near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg, name_rest
      CHARACTER*4  chari
#ifdef _MPI_
      CHARACTER*2  bchari
#endif
      CHARACTER*7  eora
      CHARACTER*40 line_res
      character*1  rshell,rshk,char_a,char_b,char_t,ato1,ato2

      INTEGER    l2,i1,j1,naus,tlscaled,atom_prot(numpr)
      INTEGER    kp,kt,hours,min,nato_slt,iter_avg,indexslice(2)
      REAL*8     wptot,xmu,xm1,xm2,tfactor
      REAL*8  uptot,upstot,ustot,ektot,E0,fact,rcom,ucom,mtar,mlig
     &     ,str_com,diffcom,r0com_t

      INTEGER ninner,nscal,ntot_fragm,fragm_1,fragm_2,ninn1,ninn0
     &     ,nsteer_temp
      INTEGER i,j,k,nsstt,nscale,ncconf,tstep,nstep,mstep,iret,ka
      INTEGER ig,flag,il,im,in0,in1,naux,nadd_tpg,iaux
      LOGICAL elflag,linit,lfirst,lprint
     &     ,lreturn,lupdate,lfalse,lrject,invertv,ladd_end,lrm_end
     &     ,laddq_end
      REAL*8  virsp_h,virs_h,virs_m,virp_h,ucns_h,ucns_l,ucns_m,ucos_h
     &     ,ucos_l,ucos_m,ucnsp_h,ucnsp_l,ucnsp_m,ucosp_h,ucosp_l
     &     ,ucosp_m,virsp_l,virs_l,virp_l,ucnp_h,ucnp_l,ucnp_m,ucop_h
     &     ,ucop_l,ucop_m,tl,tm,tn0,tn1,tm2,tm4,tl2,time2,tn12,tn02
     &     ,virsp_m,de_remove,fudgec,fstep,dnit,dip_tot,dtvi,vi(mb)
     &     ,ebussi
      REAl*8 urcs_h,urcs_l,urcs_m,urcp_h,urcp_l,urcp_m,urcsp_h,urcsp_l
     &     ,urcsp_m,virp_m,eer_h,eer_l,eer_m,stressr_h(3,3),stressr_l(3
     &     ,3),stressr_m(3,3),stressr_n1(3,3),stressr_n0(3,3)
     &     ,stressd_h(3,3),stressd_l(3,3),stressd_m(3,3),stressd_p(3,3)
     &     ,stressd(3,3),stressr(3,3),stress_conf(3,3),stress_kin(3,3)
     &     ,stress_tot(3,3),prt_m(3,3),prt_l(3,3),prt_h(3,3),prt_n1(3,3)
     &     ,prt_n0(3,3),st_m(3,3),gmgp(3,3),press_m,press_l,press_h
     &     ,press_n1,press_n0,press_conf,press_kin,velt,stresst0(3,3)

      REAL*8  temp,fcpu,elapse,puconf,pucoul,puhyd,pubnd
     &     ,ubend,uptors,uitors,uconf,ucoul,ureal,urecp,ucek,pucek,tempt
     &     ,fsrtal,fsrtalc,tempr,temppr,tcm,rcm,elaps,fnstep,rms,unbond
     &     ,cnbond,unb14,cnb14,timeq,xl,xlcut,sftalp,fsbond,fsbend
     &     ,fsin14,gsbond,gsbend,purecp,vfcp,tfcp,elps,etime,etimeq
     &     ,ttime,ttimeq,uumb,ungrp,cngrp,upconf,upcoul,timef,timefq
     &     ,gcpu,temph,uceh,hpot,ubond,aux,uslvbon,uslvben,gcpu_hd
     &     ,gcpu_ld,gcpu_md,gcpu_1,gcpu_2,gcpu_3,theoric_speed_up,ucepr
     &     ,gcpu_u,gcpu_0,gsin14,dips(3),vol_gofr,qt(4),tcmn,enisolv
     &     ,gela,gela_u,gela_hd,gela_ld,gela_md,gela_1,gela_2,gela_3
     &     ,gela_0,tela1,tela2,tcpu1,tcpu2,treal
      REAL*8  gr,pueng
      INTEGER navg,bb
      PARAMETER (navg = 60)

      REAL*8  sumarray(navg),ssmarray(navg)
      REAL*8 sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_tote,sum_ucek
     &     ,sum_temp,sum_tempt,sum_tempr,sum_temppr,sum_temph,sum_rms
     &     ,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul,sum_percip
     &     ,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor,sum_pnbd
     &     ,sum_pebnd,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     &     ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_st(3,3)
     &     ,sum_co(3,3),sum_volume,sum_pv,sum_temppra,sum_eslvint
     &     ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors

      REAL*8 ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_tote,ssm_ucek
     &     ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
     &     ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &     ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul
     &     ,ssm_epconf,ssm_presst,ssm_press,ssm_pressc,ssm_pressk
     &     ,ssm_co(3,3),ssm_st(3,3),ssm_volume,ssm_pv,ssm_temppra
     &     ,ssm_eslvint,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors
     &     ,rtime,oldtime,str_vel(500),ben_vel(500),tor_vel(500)
     &     ,bond_added(500),bend_added(500),tors_added(500),workbo(500)
     &     ,workbe(500),workto(500),Hextended,ktbeta,worktemp,wrk0,wrk0p
     &     ,dlam,dlamq,dwrkm,dwrkl,dwrkh,dwrkrc_m,dwrkrc_l,dwrkrc_h
     &     ,dwrk_al,dwrk14,dwrk14t,workcom,kcom_md

      EQUIVALENCE      
     &     (sum_econf, sumarray(1)), (sum_ecoul, sumarray(2)),
     &     (sum_enbnd, sumarray(3)), (sum_etotpot, sumarray(4)),
     &     (sum_tote, sumarray(5)), (sum_ucek, sumarray(6)),
     &     (sum_temp, sumarray(7)), (sum_tempt, sumarray(8)),
     &     (sum_tempr, sumarray(9)), (sum_temppr, sumarray(10)),
     &     (sum_temph, sumarray(11)), (sum_rms, sumarray(12)),
     &     (sum_pecek, sumarray(13)), (sum_pehyd, sumarray(14)),
     &     (sum_peconf, sumarray(15)), (sum_pecoul, sumarray(16)),
     &     (sum_enb14, sumarray(17)), (sum_ebend, sumarray(18)),
     &     (sum_ebond, sumarray(19)), (sum_eitor, sumarray(20)),
     &     (sum_eptor, sumarray(21)), (sum_pnbd, sumarray(22)),
     &     (sum_pebnd, sumarray(23)), (sum_pepot, sumarray(24)),
     &     (sum_ptote, sumarray(25)), (sum_gr, sumarray(26)),
     &     (sum_epcoul , sumarray(27)), (sum_epconf, sumarray(28)),
     &     (sum_st(1,1), sumarray(29)), (sum_co(1,1), sumarray(38)),
     &     (sum_presst, sumarray(47)), (sum_press, sumarray(48)),
     &     (sum_pressc, sumarray(49)), (sum_pressk, sumarray(50)),
     &     (sum_volume, sumarray(51)), (sum_pv, sumarray(52)),
     &     (sum_temppra,sumarray(53)), (sum_eslvint,sumarray(54)),
     &     (sum_percip,sumarray(55)),  (sum_eebond, sumarray(56)),
     &     (sum_eebend, sumarray(57)),  (sum_eeptors, sumarray(58)),
     &     (sum_eeitors, sumarray(59))
      EQUIVALENCE      
     &     (ssm_econf, ssmarray(1)), (ssm_ecoul, ssmarray(2)),
     &     (ssm_enbnd, ssmarray(3)), (ssm_etotpot, ssmarray(4)),
     &     (ssm_tote, ssmarray(5)), (ssm_ucek, ssmarray(6)),
     &     (ssm_temp, ssmarray(7)), (ssm_tempt, ssmarray(8)),
     &     (ssm_tempr, ssmarray(9)), (ssm_temppr, ssmarray(10)),
     &     (ssm_temph, ssmarray(11)), (ssm_rms, ssmarray(12)),
     &     (ssm_pecek, ssmarray(13)), (ssm_pehyd, ssmarray(14)),
     &     (ssm_peconf, ssmarray(15)), (ssm_pecoul, ssmarray(16)),
     &     (ssm_enb14, ssmarray(17)), (ssm_ebend, ssmarray(18)),
     &     (ssm_ebond, ssmarray(19)), (ssm_eitor, ssmarray(20)),
     &     (ssm_eptor, ssmarray(21)), (ssm_pnbd, ssmarray(22)),
     &     (ssm_pebnd, ssmarray(23)), (ssm_pepot, ssmarray(24)),
     &     (ssm_ptote, ssmarray(25)), (ssm_gr, ssmarray(26)),
     &     (ssm_epcoul , ssmarray(27)), (ssm_epconf, ssmarray(28)),
     &     (ssm_st(1,1), ssmarray(29)),(ssm_co(1,1),ssmarray(38)),
     &     (ssm_presst, ssmarray(47)), (ssm_press, ssmarray(48)),
     &     (ssm_pressc, ssmarray(49)), (ssm_pressk, ssmarray(50)),
     &     (ssm_volume, ssmarray(51)), (ssm_pv, ssmarray(52)),
     &     (ssm_temppra,ssmarray(53)), (ssm_eslvint,ssmarray(54)),
     &     (ssm_percip,ssmarray(55)),  (ssm_eebond, ssmarray(56)),
     &     (ssm_eebend, ssmarray(57)),  (ssm_eeptors, ssmarray(58)),
     &     (ssm_eeitors, ssmarray(59))

*=======================================================================
*    Variables needed to compute the stress tensor and implement       =
*    constant-stress simulations                                       =
*=======================================================================

      REAL*8  pressc
      REAL*8  stresstk(3,3)
      REAL*8  temppra,energy,time_fs


*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER mrject,ngrp_old,nprot_old,nind(2) ,nstep_steer,ndip
     &     ,nboxcnst,itask,istride
      INTEGER   indxi(2,indmax),indxj(2,indmax),indxk(2,indmax)

      INTEGER, allocatable :: nnlpp0(:,:),nnlpp(:,:)
     &     ,nnlpp1(:,:),cnst_protl(:,:),cnst_protlt(:,:)


#ifdef _OMP_
      integer indmax2
      parameter(indmax2=16*CELLMAX)
      INTEGER, allocatable :: nnlppf(:,:)
      INTEGER  npoints,i_ind(indmax2),j_ind(indmax2),k_ind(indmax2)
     &     ,n_ind(indmax2)
      INTEGER  nnlpp2(1)
#else
      INTEGER  nnlpp2(1),nnlppf(sitslu+1)
#endif
      REAL*8  max_dist(3),dx,dy,dz,dr
      REAL*8  fpx_n0(mb),fpy_n0(mb),fpz_n0(mb),fcax_n0(npm),fcay_n0(npm)
     &     ,fcaz_n0(npm)
      REAL*8  fpx_n1(mb),fpy_n1(mb),fpz_n1(mb),fcax_n1(npm),fcay_n1(npm)
     &     ,fcaz_n1(npm)
      REAL*8  fpx_m(mb),fpy_m(mb),fpz_m(mb),fcax_m(npm),fcay_m(npm)
     &     ,fcaz_m(npm)
      REAL*8  fpx_l(mb),fpy_l(mb),fpz_l(mb),fcax_l(npm),fcay_l(npm)
     &     ,fcaz_l(npm)
      REAL*8  fpx_h(mb),fpy_h(mb),fpz_h(mb),fcax_h(npm),fcay_h(npm)
     &     ,fcaz_h(npm)
      REAL*8  fpx(mb),fpy(mb),fpz(mb),grad_max,tipo(mb)

      REAL*8  fth(3)
      INTEGER mapdn(2,mf),nmapdn(mf)
      REAL*8  vpx(mb),vpy(mb),vpz(mb),vpx1(mb),vpy1(mb),vpz1(mb)
     &     ,vcax(numpr),vcay(numpr),vcaz(numpr),vcbx(numpr),vcby(numpr)
     &     ,vcbz(numpr),massold(mb)
      REAL*8  etap(hoov),vh1(hoov),vco(3,3),radfact

*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external

      COMMON /dynam/ fpx_h,fpy_h,fpz_h,fpx_l,fpy_l,fpz_l,fpx_m,fpy_m
     &     ,fpz_m,fpx_n0,fpy_n0,fpz_n0,fpx_n1,fpy_n1,fpz_n1,fcax_m
     &     ,fcay_m,fcaz_m,fcax_l,fcay_l,fcaz_l,fcax_h,fcay_h,fcaz_h
     &     ,fcax_n0,fcay_n0,fcaz_n0,fcax_n1,fcay_n1,fcaz_n1,vpx,vpy,vpz
     &     ,vpx1,vpy1,vpz1,etap,vh1,vcax,vcay,vcaz,vcbx,vcby,vcbz,vco
     &     ,mapdn,nmapdn

*     Phony forces and energies for neighbor list

      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p(1),fpy_p(1),fpz_p(1),rcuth_save,rtolh_save,rcutl_save
     &     ,rtoll_save,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,coul_bnd_slt_rl
     &     ,conf_bnd_slv_n1,coul_bnd_slv_n1,falch_slt
      
      LOGICAL dpress,dhoover,grflag,exist,mask(m1),h_skin,l_skin,m_skin
     &     ,n1_skin,n0_skin,time_steer,laux
      INTEGER dthree
      REAL*8 zz1(3,3),lzz,lzz1
      INTEGER igrn,krdf(maxint*g1),worka(m1,m10),type_slv(slvatm)
      INTEGER numatoms,ntype_slv,offset_slv,nbetab_slv,abmd_cryst_nvect
      INTEGER imin,jmin
      PARAMETER(nbetab_slv=90,abmd_cryst_nvect=10)
      INTEGER itype_slv(nbetab_slv),cnstco(2,5)
      CHARACTER*1 betab_slv(nbetab_slv)
      REAL*8  dssp(mbs),coeffp(mbs)
      REAL*8  xpo(mb),ypo(mb),zpo(mb),coo(3,3),co2(3,3),oc2(3,3)
      REAL*8  xpa(mb),ypa(mb),zpa(mb),xpga(m11),ypga(m11),zpga(m11)
     &     ,xpcma(npm),ypcma(npm),zpcma(npm),lx(mb),ly(mb),lz(mb)
      REAL*8  work(mspline),wca(m1),whe(m1),wbc(m1),errca(npm),errhe(npm
     &     ),errbc(npm),erral(npm),drpca(m1),drpbc(m1),drphe(m1)
     &     ,drpal(m1),xp_avg(m1),yp_avg(m1),zp_avg(m1),tmass(numpr)
     &     ,tmass1(numpr),tmassold(numpr),comscale(numpr),tmassb(numpr)
     &     ,masspp(3),dssco(5),abmd_cryst_vect(3,abmd_cryst_nvect)
     &     ,rtollo,str_tim0r,xgg,ygg,zgg,drmin,chrge0(mb)

      INTEGER cnstp(2,mbs),cnstpp,offset,abmd_dir,cnstpp_slv
     &     ,cnst_protp,count,tot_protl,tot_cnst
     &     ,nmol_slt,nmol_slv,natom_slt,natom_slv,nrestart
     &     ,n_always_accept,ierr

! saved coordinates and forces for cyclic SGE 
      real*8  xp0_s(m1),yp0_s(m1),zp0_s(m1)
      real*8  lx_s(m1),ly_s(m1),lz_s(m1)
      real*8  xpg_s(m11),ypg_s(m11),zpg_s(m11)
      real*8  xpcma_s(numpr),ypcma_s(numpr),zpcma_s(numpr)
      real*8  co_s(3,3),oc_s(3,3),vco_s(3,3),temp_s,stressr_s(3,3)
     &     ,stressd_s(3,3),gmgp_s(3,3),fth_s(3),eta_s(3)
      REAL*8  fpx_m_s(mb),fpy_m_s(mb),fpz_m_s(mb),fcax_m_s(npm)
     &     ,fcay_m_s(npm),fcaz_m_s(npm)
      REAL*8  fpx_l_s(mb),fpy_l_s(mb),fpz_l_s(mb),fcax_l_s(npm)
     &     ,fcay_l_s(npm),fcaz_l_s(npm)
      REAL*8  fpx_h_s(mb),fpy_h_s(mb),fpz_h_s(mb),fcax_h_s(npm)
     &     ,fcay_h_s(npm),fcaz_h_s(npm)
      REAL*8  fpx_n0_s(mb),fpy_n0_s(mb),fpz_n0_s(mb),fcax_n0_s(npm)
     &     ,fcay_n0_s(npm),fcaz_n0_s(npm)
      REAL*8  fpx_n1_s(mb),fpy_n1_s(mb),fpz_n1_s(mb),fcax_n1_s(npm)
     &     ,fcay_n1_s(npm),fcaz_n1_s(npm)
#ifdef _MPI_
      INTEGER ntraj_eff,para0,ibatt
#else
      INTEGER ntraj_eff,para0,ntrajectories,iproc
#endif

      COMMON /rag2/ temp,fcpu,elapse,puconf,pucoul,puhyd,pubnd,ubend
     &     ,uptors,uitors,uconf,ucoul,ureal,urecp,ucek,pucek,tempt
     &     ,fsrtal,tempr,temppr,tcm,rcm,fnstep,rms,unbond,cnbond,unb14
     &     ,cnb14,timeq,xl,xlcut,sftalp,fsbond,fsbend,fsin14,gsbond
     &     ,gsbend,purecp,vfcp,tfcp,elps,etime,etimeq,ttime,ttimeq,uumb
     &     ,ungrp,cngrp,upconf,upcoul,timef,timefq,gcpu,temph,uceh,hpot
     &     ,ubond,aux,ucepr,temppra,dssp,coeffp,xpo,ypo,zpo,xpa,ypa,zpa
     &     ,xpga,ypga,zpga,xpcma,ypcma,zpcma,wca,whe,wbc,errca,errhe
     &     ,errbc,erral,drpca,drpbc,drphe,drpal,xp_avg,yp_avg,zp_avg
     &     ,tmass,tmassb,lx,ly,lz,work,krdf,cnstp

      DATA eora/'E.O.R. '/
      DATA diffcom /0.d0/
      SAVE diffcom
      
      
*==================== EXECUTABLE STATEMENTS ============================

#ifdef _MPI_
      ntraj_eff=ntrajectories/nbatteries
#else
      ntrajectories = 1
      ntraj_eff=1
      iproc = 0
      para0 = 0
#endif

      ebussi = 0.d0
      invertv=.false.
      lfalse = .FALSE.
      IF(.NOT.clean) THEN
         write(kprint,20900)
      END IF   

*===  Check if the dimension of the work array are sufficient 

      IF(mb.LT.ntap) THEN
         errmsg='While running MD: PARAMETER MB; dim. of the work'
     &        //'arrays is insufficient. Abort. '
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF
      IF(ngrp.GT.npm) THEN  
        errmsg=' While running MD: PARAMETER NPM: dim. of the work'
     &       //' arrays is insufficient. Abort. '
        CALL xerror(errmsg,80,1,2)
        STOP
      END IF

      if (lrespa.EQ.0) lrespa=1
      if (mrespa.EQ.0) mrespa=1
      if (n1respa.Eq.0)n1respa=1
   
*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

      pi = dacos(-1.d0) 
      radfact= pi/180.d0
      nstep_steer=0

*===  set few variable to zero

      
      oldtime=0.D0
      lreturn=.false.
      lprint=.false.
      flag=nflag(1)
      IF(flag .EQ. 1) THEN
         mrject=0
         nrject=0
      END IF
      mrject=nrject*lrespa*mrespa
      time_fs=time*DFLOAT(maxstp+nrject)
      enisolv=0.d0
      nrestart=-1
      n_always_accept=0
      nscale=0
      ncconf=0
      mstep=0
      nboxcnst=5
      eer_h=0.0D0
      eer_l=0.0D0
      eer_m=0.0D0
      dip_tot=0.d0
      ndip=0
      cnbond=0.d0
      purecp=0.d0 
      coul_bnd_slt=0.d0 
      coul_bnd_slt_n1=0.d0 
      self_slt=0.d0 
      falch_slt=0.d0 
      coul_bnd_slv=0.d0 
      coul_bnd_slv_n1=0.d0 
      urcs_m=0.d0
      urcs_l=0.d0
      urcs_h=0.d0
      self_slv=0.d0
      fsrtal=0.d0 
      ninn0=0
!     initial state for OFF species (no restart) is OFF
      ladd=1.d0 
      laddq=1.d0 
!     initial state for ON species (no restart) is ON
      lrm=0.d0
      lrmq=0.d0
      ucom=0.d0
      if(addstrcom) THEN
        r0com_t=r0com
      end if
      time_steer=abs(str_time1
     &     -str_time0).gt.0.or.(path_steer.or.alchemy)
      if(time_steer.and.str_tim1r.ge.str_time1) THEN 
        str_tim0r = str_tim1r + str_time1-str_time0
      end if
      CALL zero3x3(stressr_h)
      CALL zero3x3(stressr_l)
      CALL zero3x3(stressr_m)
      CALL zero3x3(stressd_h)
      CALL zero3x3(stressd_l)
      CALL zero3x3(stressd_m)
      CALL zero3x3(stresst0)  ! this is the stress due to driven adde topology
       uumb=0.0D0
      conf_bnd_slv=0.0D0
      conf_bnd_slt=0.0D0
      tstep=0
      cnstpp=0
      fudgec=1.d0-fudge
      nsstt=0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv
      iter_avg=0
      ladd_end=.false.
      lrm_end=.false.
      laddq_end=.false.
      lrmq_end=.false.

      if(nat_added.EQ.0) ladd_end=.true.
      if(nat_removed.EQ.0) lrm_end=.true.

      IF(cpress) THEN
         DO i=1,3
            masspp(i)=masspr
         END DO
         IF(isostress) THEN 
           nboxcnst=5
           CALL set_const_co(co,dssco,cnstco)
         ELSE IF(isostressxy) THEN  
           nboxcnst=4
           CALL set_constxy_co(co,dssco,cnstco)
           isostress=.true.
         END if
      END IF
      CALL comp_molmass(nprot,protl,mass,tmass)

      IF(abmd) THEN
         IF(fold) THEN
            IF(.NOT. abmd_unbias) THEN
               IF(rspset .LT. 0.0D0) THEN
                  abmd_dir=-1
               ELSE
                  abmd_dir=1
               END IF
               rspset=0.0D0
            END IF
         ELSE IF(dissociate .OR. associate) THEN
            IF(dissociate) THEN
               abmd_dir=1
            ELSE IF(associate) THEN
               abmd_dir=-1
            END IF
            rspset=0.0D0
            CALL find_diss_mol(diss_mol,diss_atoms,atoms_diss,mol_diss
     &           ,diss_list,nprot,protl)
         ELSE IF(abmd_cryst) THEN
            rspset=0.0D0
            abmd_dir=abmd_cryst_dir
            CALL gen_abmd_kvect(abmd_kvect,abmd_cryst_nvect
     &           ,abmd_cryst_vect)
         END IF
      END IF

      if(meta_run) then
        if(metar) then      
! do not read anything but the number of old hills and positions: !!!!
! real parameters come from the input file
          read(kmetar,*) meta_dc, meta_n
          do i = 1,meta_n
            read(kmetar,*) (meta_list(i,j), j = 1,meta_nrc)
          enddo
        endif
! open metadynamics output
        meta_data_file='META_DATA'
        call openf(kmeta,meta_data_file,'FORMATTED','UNKNOWN',0)
! convert the height of hills to program units 
        meta_height = meta_height * 1000. / (unite*avogad)
      endif
      
*=======================================================================
*--- Initialize the velocities if at the beginning of the run        ---
*=======================================================================

      CALL zeroa(vpx,vpy,vpz,ntap,1)
      CALL zeroa(vpx1,vpy1,vpz1,ntap,1)
      CALL zeroa(vcax,vcay,vcaz,nprot,1)
      CALL zeroa(xp_avg,yp_avg,zp_avg,ntap,1)

*=======================================================================
*----- Initialize field arrays if needed 
*=======================================================================

      if(ncofactor.gt.0.and.nvi.gt.0) then
c--      set field cut-off to outer shell radius and alphaf 
c--      to current value if non set    
         if(alphaf.lt.0.d0.or.cut_field.lt.0.d0) then
            alphaf=alphal
            cut_field=rcuth+rtolh
         end if
         if(ncofactor.ne.12345678) THEN 
           CALL field_init(ncofactor,nsegs,kprint,naux) 
           if(naux.eq.0) THEN
             ncofactor = 0
             write(kprint,10977) 
           endif
         ELSE
           ncofactor=ntap
           do i=1,ntap 
             nsegs(1,i)=i
             nsegs(2,i)=i
           end do
         ENDIF
         dtvi = time/float(lrespa*mrespa)
         cut_field=cut_field**2
         IF(flag.eq.0.or.flag.eq.2) THEN  
           if(lform_field) THEN 
             write(kvi,10103) dtvi,ncofactor
10103        FORMAT(f20.5,I10)
             if(ncofactor.ne.ntap) THEN 
               do i=1,ncofactor
                 write(kvi,10104) nsegs(1,i),nsegs(2,i)
10104            FORMAT(2I10) 
               end do
             END IF
           ELSE
             write(kvi) dtvi,ncofactor
             do i=1,ncofactor
               write(kvi) nsegs(1,i),nsegs(2,i)
             end do
          END IF
         END IF
      end if

      if(start_conf) then 
         max_dist(1) = distmax
         max_dist(2) = distmax
         max_dist(3) = distmax
      end if

      IF(nplot_fragm .GT. 0) THEN
         do i=1,nfragm
            IF ( fragm(1,i).gt.ntap.or.fragm(2,i).gt.ntap) THEN 
               write(kprint,10260) i
10260          FORMAT(//  
     &              '**ERROR: fragments limits for',i5,' DEF_FRAGM',
     &              '(&READ_SOLUTE) exceeds protein atoms' )
               STOP
            END IF 
         end do
         ntot_fragm=0
         do i=1,nfragm
            fragm_1=fragm(1,i)
            fragm_2=fragm(2,i)
            ntot_fragm = ntot_fragm + fragm_2-fragm_1+1 
         end do    
      END IF

*=======================================================================
*----- allocate neighbor lists
*=======================================================================

#ifdef _OMP_
      mpp8=1+mpp/nthr
      allocate(nnlppf(mpp8,nthr),stat=iret)
      allocate(nnlpp0(mpp8,nthr),stat=iret)
      allocate(nnlpp(mpp8,nthr),stat=iret)
      allocate(nnlpp1(mpp8,nthr),stat=iret)
      m1t=1+2*m1/nthr2
      allocate(cnst_protl(2*m1,1),stat=iret)
      allocate(cnst_protlt(m1t,nthr2),stat=iret)
#else
      allocate(nnlpp0(mpp,1),stat=iret)
      allocate(nnlpp(mpp,1),stat=iret)
      allocate(nnlpp1(mpp,1),stat=iret)
      allocate(cnst_protl(2*m1,1),stat=iret)
#endif


*=======================================================================
*----- Copy old CO and OC matrix to temporary arrays -------------------
*=======================================================================

      IF(change_cell) THEN
         CALL dcopy(9,co,1,co2,1)
         CALL dcopy(9,oc,1,oc2,1)
      END IF

*=======================================================================
*----- Read input from restart file if needed --------------------------
*=======================================================================

      kp=0
      kt=-1
      IF(flag.GT.0) THEN
         elflag=electr.AND.(.NOT.elinit)
         write(kprint,70300) 
         CALL readrs(kdump,nstep,temp,ntap,ngrp_old,nprot_old,xp0,yp0
     &        ,zp0,vpx,vpy,vpz,xpg,ypg,zpg,xpcma,ypcma,zpcma,vpx1,vpy1
     &        ,vpz1,vcax,vcay,vcaz,eta,etap,vh1,dthree,dhoover,sumarray
     &        ,ssmarray,navg,navg,co,oc,dpress,vco,zz1,lzz,lzz1,grflag
     &        ,krdf,igrn,maxint*g1,restart_cont)
         IF(hmass) CALL chmass(beta,mass,ntap,hdmass)
         write(kprint,*) "     -------"
     &        ,TRIM(restart_file ),"--------"
         write(kprint,70400)
#ifdef _MPI_
!        re-inizialize the velocities if battery-restart.
!        Each battery will follow and independent GE trajectory.         
         if(nbatteries.gt.1) THEN
           write(kprint,10081) 
10081      FORMAT(//,
     &         5x, '<----  REM run with batteries -------------->',/
     &         5x,'<      Velocities re-inizialized            >',/ 
     &         5x,'<-------------------------------------------> ')
           LINIT=.true.
           CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,linit
     &          ,massinfty)
           CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &          ,vcaz)
           LINIT=.FALSE.
           IF(cpress) THEN
             CALL set_tempp(masspr,vco,temppra,0.0D0)
           END IF
           IF(thermos) THEN
             CALL set_tempt(neta,qmass,etap,temph,0.0D0)
           END IF
         end if
#endif         
10091    ninner=nstep*lrespa*mrespa
         ninn1=nstep*lrespa*mrespa*n1respa

c-- check whether one atom has mass > 1D20 and normal velocities.
c-- this can happen restarting a SMD using tpgprm with large masses and a
c-- restart file produced with a normal masses.   
         do i=1,ntap 
           if(mass(i)*abs(vpx(i)).GT.1.D10.and.mass(i)
     &          *abs(vpy(i)).GT.1.D10.and.mass(i)*abs(vpz(i)).GT.1.D10)
     &          THEN
             vpx(i)=0.d0 
             vpy(i)=0.d0 
             vpz(i)=0.d0
             write(kprint,17543) i,betb(i)
17543        FORMAT("   ---- Atom ",i5,2x,a7
     &            ," is fixed: velocities set to 0")

           ENDIF
         end do
c--      RESTART A SMD simulation 
         if(flag.eq.1.and.time_steer) THEN
c           go to 9998
           if(path_steer) THEN 
             rtime=nstep*time
c--          this is done for restarting a steer_path SMD.
             if(.not.alchemy) THEN 
               call steer_along_path(lstretch,nbonds_added,lbend
     &              ,nbends_added,litor,nitors_added,strbonds
     &              ,strbends,strtors,potbo(1,2),potbe(1,2)
     &              ,potit(1,2),potbo(1,1),potbe(1,1)
     &              ,potit(1,1),radfact,timerec,ntimes,pathbo
     &              ,pathbe,pathto,activebo,activebe,activeto,rtime
     &              ,eqdist,eqang,eqdied,str_vel,ben_vel,tor_vel
     &              ,indexslice)
               ELSE
                 do i=1,ntap
                   ss_alchemy(i)=1
                 end do
                 call initialize_lambda(nat_added,nat_removed,ladded
     &                ,laddedq,lremoved,lremovedq,atom_added
     &                ,atom_removed,lambda,lambdaq,listqq,nlistqq
     &                ,lstretch,lstrtch,lbend,lbndg,int14p,int14,lconstr
     &                ,lcnstr,clewld,ss_alchemy)
                 call alchemic(timerec,ntimes,ladded,laddedq,ladd,laddq
     &                ,lremoved,lremovedq,lrm,lrmq,rtime,indexslice)
                 ladd0=ladd
                 laddq0=laddq
                 lrm0=lrm
                 lrmq0=lrmq
                 if(ladd.NE.0.D0.OR.laddq.NE.0.d0) THEN 
                   do i=1,nat_added
                     lambda(atom_added(i))=ladd
                     lambdaq(atom_added(i))=laddq
                     lambda0(atom_added(i))=ladd
                     lambdaq0(atom_added(i))=laddq
                   end do
                 ENDIF
                 if(lrm.NE.0.D0.OR.lrmq.ne.0.d0) THEN 
                   do i=1,nat_removed
                     lambda(atom_removed(i))=-lrm
                     lambdaq(atom_removed(i))=-lrmq
                     lambda0(atom_removed(i))=-lrm
                     lambdaq0(atom_removed(i))=-lrmq
                   end do
                 ENDIF
               END IF 
           ELSE
c--          this is done for restarting a standard SMD
             do i=1,nbonds_added
               if(strbonds(i).EQ.1) THEN
                 aux = (eqdist1(i)-eqdist(i)) /(str_time1 - str_time0)
c-               recalculate initial values of eqdist
                 eqdist(i)=eqdist(i)+ aux*(nstep*time-str_time0) 
               END IF
             end do
             do i=1,nbends_added
               if(strbends(i).EQ.1) THEN
                 aux = (eqang1(i)-eqang(i)) /(str_time1 - str_time0)
c-               recalculate initial values of eqang
                 eqang(i)=eqang(i)+ aux*(nstep*time-str_time0) 
               END IF
             end do
             do i=1,nitors_added
               if(strtors(i).EQ.1) THEN 
                 aux = (eqdied1(i)-eqdied(i)) /(str_time1 - str_time0)
c-               recalculate initial values of eqdied
                 eqdied(i)=eqdied(i)+ aux*(nstep*time-str_time0) 
               END IF
             end do
             str_time0=nstep*time
             if(addstrcom) THEN 
               aux = (r1com-r0com) /(str_time1 - str_time0)
c-             recalculate initial values of r0com
               r0com_t=r0com+aux*(nstep*time-str_time0) 
             end if
           end if
9998       continue
         end if

!   compute restart time 
         IF(flag .EQ. 1) THEN
            time_fs=time*DFLOAT(maxstp-nstep)
            oldtime=nstep*time
         END IF
      END IF

c--   Scale masses of atoms if command SCALE_MASS (&RUN) is given
      if(scalemass) THEN 
        CALL scale_masses(ntap,mass,atomtoscale,velscale,nscalemass,vpx
     &       ,vpy,vpz)
        CALL comp_molmass(nprot,protl,mass,tmass)
c--        Computed new com velocities
        CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &       ,vcaz)
#ifdef _OMP_
        call assign_prot(nprot,protl,ss_index,iprot,typei,tmass
     &       ,tmass1)
#endif
        write(kprint,10988) 
10988   FORMAT(//'==== MASSES HAVE BEEN SCALED ================' //
     &       15x,' atom',15x,' scalefact')
        do i=1,nscalemass 
          k=atomtoscale(i)
          write(kprint,10999) k,beta(k),1/velscale(i)**2
10999     FORMAT(5x,i10,2x,A7, "  ----->  ", 2G15.5)  
          aux=mass(k)/velscale(i)**2.
          if(aux.gt.massinfty) massinfty=2*aux
        end do
      END IF

      IF(flag .EQ. 1 .AND. nstep .GE. maxstp ) THEN
         WRITE(errmsg,'(''In MTSMD: Step number in restart = '',
     &        i6,'' while in input = '',i6,'' Increase TIME.'')')
     &        nstep,maxstp
         CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*---  Write a Banner Page                                            ---
*=======================================================================

      write(kprint,'(//)')
      write(kprint,1200)
      write(kprint,1300)
      write(kprint,1100)
      write(kprint,1100)
      write(kprint,1300)
      write(kprint,1200)
      write(kprint,'(//)')

*=======================================================================
*---  Reset the averages when NFLAG = 0 or 2                         ---
*=======================================================================

      IF(flag.EQ.2. .OR. flag .EQ. 3 .OR.flag.EQ.0 
     &     .OR. flag .EQ. 4) THEN
          CALL zero(sumarray,ssmarray,navg,navg)
         nstep=0
         ninner=0
         ninn1=0
         ninn0=0
      END IF

      IF(anxrms .OR. gofr .OR. avg_str) CALL asng_xrmsw(ss_point,m1+1
     &     ,wca,whe,wbc,beta,mback,nbone)

      IF(anxrms_cell) THEN
         anprot=.TRUE.
         annpro=annpro+1
         anpoint(1,annpro)=1
         anpoint(2,annpro)=ntap
      END IF


      IF(gofr) CALL get_type_slv(nato_slt,nato_slv,beta,betab_slv
     &     ,nbetab_slv,ntype_slv,type_slv,itype_slv,offset_slv
     &     ,types_gofr,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*---- Initialize REM simulation ----------------------------------------
*=======================================================================
#ifdef _MPI_
      if( rem_run ) call init_rem(iproc,ntrajectories,kprint,ntap
     &     ,ltors,litor,int14p,ltor,litr,int14,remltor,remitor,remint14
     &     ,rematom)
#endif
*=======================================================================
*---- Initialize SGE simulation ----------------------------------------
*=======================================================================
      if( sim_tempering ) call sge_init( flag, iproc, ntrajectories,
     &     nbonds_added, nbends_added, nitors_added, ntap, ltors,
     &     int14p, ltor, int14, xp0, yp0, zp0, radfact, force_cost,
     &     eqdist, eqdist1, force_ang, eqang, eqang1, force_died,
     &     eqdied, eqdied1, atom_b1, atom_b2, atom_be1, atom_be2,
     &     atom_be3, atom_it1, atom_it2, atom_it3, atom_it4, strbonds,
     &     strbends, strtors, remltor, remint14, rematom )

*=======================================================================
*---- Initialize G of Rs -----------------------------------------------
*=======================================================================

      IF(gofr) THEN
         CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
      END IF

*=======================================================================
*-------- Compute two vectors used by the constraint subroutine --------
*=======================================================================

      write(kprint,'(//,80("="),/,20x,"FORCE FIELD IS --> ",A4,/,80("
     &=")//   )')ffield_type(1:4)

      IF(lconstr .NE. 0) THEN
         IF(adjust_cnstr) THEN
            write(kprint,14000) 
            CALL adjust_bonds(ss_index,ntap,lcnstr,lconstr
     &        ,xp0,yp0,zp0,potbo_cnst,M9,iret,errmsg,ndim)
            IF(iret .EQ. 2) CALL xerror(errmsg,80,1,20)
         END IF
         cnstpp=0
         CALL bcnstp(ss_index,lcnstr,lconstr,xp0,yp0,zp0,mass,
     &        dssp,coeffp,cnstp,mbs,cnstpp,cnstpp_slv,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         IF(stretch_heavy) THEN
            CALL get_prot_cnstr(ntap,clsthl,nclsth,cnstp,cnstpp
     &           ,cnst_protl,cnst_protp,tot_protl,tot_cnst,mask,worka
     &           ,m1)
#ifdef _OMP_
      if(omp_dynamic) call omp_set_dynamic(.TRUE.)

!$OMP PARALLEL SHARED(istride) PRIVATE(itask) NUM_THREADS(nthr2)
            itask=1+OMP_GET_THREAD_NUM()
            do i=1,m1t
              cnst_protlt(i,itask)=0
            end do
!           Build constraint map for thread ITASK
            call split_cnstr(cnst_protp,cnst_protl,cnst_protlt(1
     &           ,itask),cnstp_threads,itask,nthr2)
!$OMP END PARALLEL 
            write(kprint,'(//)')
            if(omp_dynamic) THEN 
              call omp_set_dynamic(.FALSE.)
              call OMP_SET_NUM_THREADS(nthr) !restore thread numebr
            END IF
            deallocate(cnst_protl,stat=iret) ! release memory for the big map 
            if ( iret /= 0 ) stop 'Deallocation failed for cnst_protl'
            allocate(cnst_protl(m1t,nthr2)) ! re-allocate cnst map array 
            do i=1,nthr2
              do j=1,cnstp_threads(i)
                cnst_protl(j,i)=cnst_protlt(j,i) ! copy temp map on cnst map
              end do
            end do
            deallocate(cnst_protlt,stat=iret) ! release memory  for tmp map
            if ( iret /= 0 ) stop 'Deallocation failed for cnst_protlt'
#endif            
         ELSE
           if(nthr.gt.1) THEN 
             errmsg='Use STRETCHING HEAVY with OpenMP.'
     &     //'. If solute is large, load-balance is poor'
             CALL xerror(errmsg,80,1,1)
             STOP
           ENDIF
           CALL get_prot_cnstr(ntap,protl,nprot,cnstp,cnstpp
     &          ,cnst_protl,cnst_protp,tot_protl,tot_cnst,mask,worka
     &            ,m1)
         END IF
      ELSE
        cnstpp_slv=0
      END IF


      CALL comp_nmol(ss_index,nprot,protl,nmol_slt,nmol_slv,natom_slt
     &        ,natom_slv)

#ifdef _OMP_
!     compute vectors for xp0 recalculation from lx0 in parallel
      call assign_prot(nprot,protl,ss_index,iprot,typei,tmass,tmass1)
#endif

*================================================================================
!---  setup array of fixed atoms (mass > 1.d30) 
*================================================================================

      j=0
      nfixed=0
      do i=1,ntap
        if(mass(i).gt.massinfty) THEN 
          j=j+1
          ifixed(j)=i
        end if  
      end do
      nfixed=j

*================================================================================
!---  Move dangerous (fast) proper torsion to improper  
*================================================================================
      
      tfactor=4.d0*4.184d3*1.d20/0.001
      if(moveitors) THEN
        do i=1,ltors
          xm1=12.0
          xm2=12.0
          ato1 = betb(ltor(1,i))
          if(ato1.eq."*") ato1=betb(ltor(1,i))(2:2) 
          ato2 = betb(ltor(4,i))
          if(ato2.eq."*") ato2=betb(ltor(4,i))(2:2) 
          if(ato1.eq."h") xm1=1.
          if(ato2.eq."h") xm2=1.
          wptot=abs(potto(i,1))*unite*avogad/(4.184*1000.d0) ! back to kcal  
          naux = nint(potto(i,2)) !multiplicity 
!       compute the torsional uncoupled frequency with 
!       d set to 1 (worse case scenario) 
          xmu=(xm1+xm2)*sin(3.14d0/3.d0)**2 
          aux =((tfactor*wptot/xmu)**0.5)*dfloat(naux)/(3.d10*2*3.14)
          if(aux.gt.300.d0) THEN   
            litor=litor+1  
            potit(litor,1)=potto(i,1)
            litr(1,litor)=ltor(1,i)
            litr(2,litor)=ltor(2,i)
            litr(3,litor)=ltor(3,i)
            litr(4,litor)=ltor(4,i)
            potit(litor,2)=potto(i,2)
            potit(litor,3)= -1.d0 ! cosine improper torsion
            potto(i,1) = 0.d0   ! torsion is set to zero and will be skipped in fptors
            if(rem_run.AND.rem_segment) THEN ! check if torsions are in rem_segment if rem_run  
              do j = 1,rem_groups
                do i1 = 1,4
!               if the atoms of the torsion overlap with the subsets of
!               atoms, then 
                  if((litr(i1,litor).GE.rem_group(1,j)).AND.(litr(i1
     &                 ,litor).LE.rem_group(2,j))) remitor(litor) =
     &                 .TRUE.
                enddo
              enddo
            end if
          end if
        end do  
      end if
*================================================================================
!---  add extra bonds/bend/tors  given in READ_POTENTIAL using ADD_STR_BONDS cmd
*================================================================================
      
      write(kprint,20018)
20018 FORMAT(20("=")," Added Topology (driven or not) starts here",
     &     17("="))   
      if(nbonds_added.gt.0) THEN 
        write(kprint,1210) nbonds_added
1210    FORMAT (/ ' -----  Added ',i4, ' extra bonds (MTSMD LEVEL) '//
     &       13x,4x,'atom1',4x,'atom2',5x,'K',8x,'r0',8x,'r1',8x
     &       ,'vel(A/ns)')
        write(kprint,20008) 
        do i=1,nbonds_added
          if(atom_b1(i).le.nato_slt) THEN
            lstrtch(1,i+lstretch)=atom_b1(i)
          ELSE
            errmsg = 'ADD_BONDS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_b2(i).le.nato_slt) THEN
            lstrtch(2,i+lstretch)=atom_b2(i)
          ELSE
            errmsg = 'ADD_BONDS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          ENDIF
          if(strbonds(i).EQ.1.and.time_steer) THEN
            workbo(i)=0.d0
            potbo(i+lstretch,1) =  1000.0D0*force_cost(i)*4.184/(unite
     &           *avogad)
            str_vel(i) = (eqdist1(i)-eqdist(i)) / (str_time1 -
     &           str_time0)
            potbo(i+lstretch,2) =  eqdist(i)
            write(kprint,1211) lstrtch(1,i+lstretch)
     &           ,lstrtch(2,i+lstretch),force_cost(i),eqdist(i)
     &           ,eqdist1(i),1000000.d0*str_vel(i)

1211        FORMAT(10x,2i10,4g10.3)
          ELSE
            potbo(i+lstretch,1) =  1000.0D0*force_cost(i)*4.184/(unite
     &           *avogad)
            potbo(i+lstretch,2) =  eqdist(i)
            write(kprint,1201) lstrtch(1,i+lstretch)
     &           ,lstrtch(2,i+lstretch),force_cost(i),eqdist(i)
1201        FORMAT(10x,2i10,2f10.3,2(9x,'-'))
          END IF

        end do
        lstretch=lstretch+nbonds_added
      END IF

c---  add extra bends given in READ_POTENTIAL using ADD_BENDINGS cmd
      
      if(nbends_added.gt.0) THEN 
        write(kprint,1310) nbends_added
1310    FORMAT (/ ' --  Added ',i4, ' extra bendings (MTSMD LEVEL) '/
     &       13x,4x,'atom1',4x,'atom2',4x,'atom3',5x,'K',8x,'A0',8x
     &       ,'A1',8x,'v(deg/ns)')
        do i=1,nbends_added
          if(atom_be1(i).le.nato_slt) THEN
            lbndg(1,i+lbend)=atom_be1(i)
          ELSE
            errmsg = 'ADD_BENDINGS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_be2(i).le.nato_slt) THEN
            lbndg(2,i+lbend)=atom_be2(i)
          ELSE
            errmsg = 'ADD_BENDING atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          ENDIF
          if(atom_be3(i).le.nato_slt) THEN
            lbndg(3,i+lbend)=atom_be3(i)
          ELSE
            errmsg = 'ADD_BENDINGS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          ENDIF
          if(strbends(i).EQ.1.and.time_steer) THEN
            workbe(i)=0.d0
            potbe(i+lbend,1) = 1000.0D0*force_ang(i)*4.184/(unite
     &           *avogad)
            ben_vel(i) = radfact*(eqang1(i)-eqang(i)) / (str_time1 -
     &           str_time0)
            potbe(i+lbend,2) =  eqang(i)*radfact
            write(kprint,1301) lbndg(1,i+lbend),lbndg(2
     &           ,i+lbend),lbndg(3,i+lbend),force_ang(i),eqang(i)
     &           ,eqang1(i),1000000.d0*ben_vel(i)/radfact
1301        FORMAT(10x,3i10,4g10.3)
          ELSE
            potbe(i+lbend,1) = 1000.0D0*force_ang(i)*4.184/(unite
     &           *avogad)
            potbe(i+lbend,2) =  eqang(i)*radfact
            write(kprint,1302) lbndg(1,i+lbend),lbndg(2
     &           ,i+lbend),lbndg(3,i+lbend),force_ang(i),eqang(i)
1302        FORMAT(10x,3i10,2f10.3,2(9x,'-'))
          END IF
        end do
        lbend=lbend+nbends_added
      end if

c---  add extra itors given in READ_POTENTIAL using ADD_TORSIONS cmd
      
      if(nitors_added.gt.0) THEN 
        write(kprint,1410) nitors_added
1410    FORMAT (/ ' --  Added ',i4, ' extra torsions (MTSMD LEVEL) '/
     &       13x,4x,'atom1',4x,'atom2',4x,'atom3',4x,'atom4',5x,'K',8x
     &       ,'A0',8x,'A1',8x,'v(deg/ns)')
        do i=1,nitors_added
          if(atom_it1(i).le.nato_slt) THEN
            litr(1,i+litor)=atom_it1(i)
          ELSE
            errmsg = 'ADD_TORSIONS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_it2(i).le.nato_slt) THEN
            litr(2,i+litor)=atom_it2(i)
          ELSE
            errmsg = 'ADD_TORSIONS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_it3(i).le.nato_slt) THEN
            litr(3,i+litor)=atom_it3(i)
          ELSE
            errmsg = 'ADD_TORSIONS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_it4(i).le.nato_slt) THEN
            litr(4,i+litor)=atom_it4(i)
          ELSE
            errmsg = 'ADD_TORSIONS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          IF(strtors(i).EQ.1.and.time_steer) THEN 
            workto(i)=0.d0
            potit(i+litor,2)=eqdied(i)*radfact
            potit(i+litor,1)=1000.0D0*force_died(i)*4.184/(unite
     &           *avogad)
            potit(i+litor,3) = 1.d0
            tor_vel(i) = radfact*(eqdied1(i)-eqdied(i)) / (str_time1 -
     &           str_time0)
            write(kprint,1401) litr(1,i+litor),litr(2,i
     &           +litor) ,litr(3,i+litor),litr(4,i+litor),force_died(i)
     &           ,eqdied(i),eqdied1(i),tor_vel(i)*1000000.d0/radfact
1401        FORMAT(10x,4i10,4g10.3)
          ELSE
            potit(i+litor,2)=eqdied(i)*radfact
            potit(i+litor,1)=1000.0D0*force_died(i)*4.184/(unite
     &           *avogad)
            potit(i+litor,3) = 1.d0
            write(kprint,1402) litr(1,i+litor),litr(2,i
     &           +litor) ,litr(3,i+litor),litr(4,i+litor),force_died(i)
     &           ,eqdied(i)
1402        FORMAT(10x,4i10,2f10.3,2(9x,'-'))
          end if
        end do
        write(kprint,20008) 
        litor=litor+nitors_added
      end if
      if(addstrcom) THEN 
        write(kprint,1422) ilig1,ilig2,itar1,itar2,kcom,r0com
1422    FORMAT (/ ' -- Added COM restraint potential ---------' //
     &       10x,"Ligand: ",i10, ' ----> ', i10 / 
     &       10x,"Target: ",i10, ' ----> ', i10 / 
     &       10x,"kcom  ", f12.5, ' kcal mol-1 Angs-2' /
     &       10x,"r0com ", f12.5, ' Angs' /) 
        kcom_md=1000.0D0*kcom*4.184/efact ! convert to MD units
        if(steer_com) THEN
          workcom=0.d0
          str_com = (r1com-r0com) / (str_time1 -
     &           str_time0)
          write(kprint,1423) str_time0,str_time1,1000000.d0*str_com
1423      FORMAT(10x, " COM Steering is ON " /
     &           10x, "Start Time", G12.5, /  
     &           10x, "Final Time", G12.5, /
     &           10x, "steering velocity ", G12.5, ' Angs/ns ')
        end if
        write(kprint,20019)
      endif 

      if(path_steer) THEN 
        if(.not.alchemy) THEN
          write(kprint,20002) nbonds_steer_added
     &         ,nbends_steer_added,nitors_steer_added 
20002     format(/79("-"),/,20(" "), "STEER PATH ACTIVE with: ",/,
     &         25(" "), "steered nbonds =",i2, /,
     &         25(" "), "steered nbends =",i2, /,
     &         25(" "), "steered tors   =",i2,/)
          write(kprint,20005) 
20005     FORMAT(20x,4("="),"PATH DEFINITION AND TIME PROTOCOL",4("="),
     &         /) 
          write(kprint,20003) 
20003     FORMAT("Time Slice  Time/ps  ====== Driven Coordinates (bonds(
     &B),""bends(A),tors(T))"/79("-"))
          
          char_b="B"
          char_a="A"
          char_t="T"
          write(kprint,20007) (char_b,i=1
     &         ,nbonds_steer_added),(char_A,i=1,nbends_steer_added)
     &         ,(char_t,i=1,nitors_steer_added)  
20007     FORMAT(20x,15(4x,A1,4x))
          do i=1,ntimes
            write(kprint,20004)  i,timerec(i)/1000.d0,
     &           (pathbo(i,j),activebo(i,j),j=1,nbonds_steer_added),
     &           (pathbe(i,j),activebe(i,j),j=1,nbends_steer_added),
     &           (pathto(i,j),activeto(i,j),j=1,nitors_steer_added)
20004       FORMAT(5x,i3,2x,f7.1,2x,15(f7.2,i2))
          end do
          write(kprint,20008) 
20008     FORMAT(79("-"))
          write(kprint,20019)
20019     FORMAT(20("=")," Added Topology ends here ",33("="))   
        ELSE
!         check whether alchemical work has been read in from the
!         restart file
          if(.not.NEAR0(abs(wrk_alchemy))) wrk_alchemy=wrk_alchemy
     &         *1000.d0/efact
          if(kplot_alch0.gt.0) write(kplot_alch0,40002)
40002     FORMAT("#",7x,"Time/ps",6x,"lamb-added",13x,"lamb-rmvd",8x
     &         ,"    work", / 
     &           '#',6x,13x,"LJ",8x,"QQ",11x,"LJ",8x,"QQ")
          write(kprint,30003) nat_added,nat_removed
     &         ,do_workq,print_qwork,intramol,ljintra 
30003     format(/79("-"),/,20(" "),  
     &         "====ALCHEMICAL TRANFORMATION===== ",/,
     &         25(" "), "atom to be switched on   =",i4, /,
     &         25(" "), "atom to be switched off  =",i4, /,
     &         25(" "), "Ewald work               =",L2, /,
     &         25(" "), "Print Ewald work         =",L2  /,
     &         25(" "), "Intramol                 =",L2  /,
     &         25(" "), "LJ on                    =",L2  )
!         set lambda shift values and check whether alchemical definition is OK 
!         set initial values of lambdas for cold restart
          if(flag.ne.1) THEN
            if(nat_added.gt.0) THEN 
              linitq= .false. 
            END IF
            do i=1,ntap
              lambda(i)=0.d0 
              lambdaq(i)=0.d0 
            end do
            do i=1,ntap
              ss_alchemy(i)=1
            end do
            call initialize_lambda(nat_added,nat_removed,ladded
     &                ,laddedq,lremoved,lremovedq,atom_added
     &           ,atom_removed,lambda,lambdaq,listqq,nlistqq,lstretch
     &           ,lstrtch,lbend,lbndg,int14p,int14,lconstr,lcnstr
     &           ,clewld,ss_alchemy)
            do i=1,ntap
              lambda0(i)=lambda(i)
              lambdaq0(i)=lambdaq(i)
            end do
         END IF
          
!         initialize shift values for atoms 
          do i=1,ntap
            rmin(i)=pnbd1(nbtype(i))/2.d0**(1./6.)
!           abs below is to avoid problems with minus sign in mixing rule
            epsm(i)=abs(pnbd2(nbtype(i))*(1000.d0*4.184/efact)) 
            if(epsm(i).eq.0.d0) rmin(i)=0.1
          end do  
          do i=1,nat_added
            if(atom_added(i).gt.nato_slt) THEN 
              errmsg =
     &        'DEFINE_ALCHEMICAL ATOMS out of solute range (added)'
              CALL xerror(errmsg,80,1,1)
              STOP
            END IF
          end do
          do i=1,nat_removed
            if(atom_removed(i).gt.nato_slt) THEN 
              errmsg = 
     &        'DEFINE_ALCHEMICAL ATOMS out of solute range (removed)'
              CALL xerror(errmsg,80,1,1)
              STOP
            END IF
          end do
          write(kprint,20005) 
          write(kprint,30004) 
30004     FORMAT("Time Slice   Time/ps  ====== l(add)LG  l(add)QQ"
     &         ,"  l(Rm)LG l(Rm)QQ  "/79("-"))
          
          do i=1,ntimes
            write(kprint,30005)  i,timerec(i)/1000.d0
     &           ,ladded(i),laddedq(i),lremoved(i),lremovedq(i)
30005       FORMAT(5x,i3,2x,f8.2,11x,2(2f8.3,2x))
          end do
          write(kprint,20008) 
          IF(shifted) THEN
            write(kprint,'(10x,"SHIFT-POTENTIALS -> Sh_elec"," =",f8.3
     &,"    Sh_LJ =",f8.3 )') qshift,ljshift
          ELSE
            write(kprint,'(10x,"SOFT-CORE  -> SC_elec"," =",f8.3
     &,"    SC_LJ =",f8.3 )') qsoft,ljsoft
          ENDIF  
          write(kprint,20019)
        END IF
      end if 


*================================================================================
c---  end of adding extra bonds/bend/tors  
*================================================================================


*================================================================================
c--   Check whether initial values for steered parameters are compatible
c--   with initial actual values of the associated driven coordinate 
*================================================================================

      if(addstrcom) THEN 
        laux=.false.
        aux=0.d0
        do i=ilig1,ilig2
          aux=aux+mass(i)
        end do
        mlig=1.d0/aux
        aux=0.d0
        do i=itar1,itar2
          aux=aux+mass(i)
        end do
        mtar=1.d0/aux

        call com_restraint (kcom_md,r0com,xp0,yp0,zp0,mlig,mtar,mass
     &       ,rcom,fpx_p,fpy_p,fpz_p,ucom,ilig1,ilig2,itar1,itar2,laux)
        if(abs(kcom*(r0com-rcom)).gt.40.d0) THEN 
          write(kprint,7564)r0com,rcom 
7564      FORMAT("****** WARNING ! WARNING ! WARNING: " //
     &         " starting value for driven of COM distance is ",
     &         f10.4, " while actual value is ", f10.4 //    
     &         "****** WARNING ! WARNING ! WARNING: " //)    
        ELSE
          write(kprint,7565) r0com,rcom 
7565      FORMAT(//
     &         " starting value for driven of COM distance is ",
     &         f10.4, ".   Actual value is ", f10.4 // )   

        END IF
      END IF

      nadd_tpg=nbonds_added+nbends_added+nitors_added
      if(nadd_tpg.gt.0) call steer_variables(nbonds_added,nbends_added
     &     ,nitors_added,atom_b1,atom_b2,atom_be1,atom_be2,atom_be3
     &     ,atom_it1,atom_it2,atom_it3,atom_it4,xp0,yp0,zp0,bond_added
     &     ,bend_added,tors_added,potit(1,2),litor,nitors_added)


      if(nbonds_added.gt.0) THEN 
        do i=1,nbonds_added
          if(abs(force_cost(i)*(eqdist(i)-bond_added(i))).gt.40.d0) THEN 
            write(kprint,7544) i,eqdist(i),bond_added(i) 
7544        FORMAT("****** WARNING ! WARNING ! WARNING: " //
     &           " starting value for driven", i4,"-th stretching is ",
     &           f10.4, " while actual value is ", f10.4 //    
     &           "****** WARNING ! WARNING ! WARNING: " //)    
          END IF
        end do
      end if 
      if(nbends_added.gt.0) THEN 
        do i=1,nbends_added
          aux = abs(force_ang(i)*(eqang(i)-bend_added(i)*180./pi))
          if(aux.gt.400.d0) THEN 
            write(kprint,7545) i,eqang(i),bend_added(i)
     &           *180.d0/pi 
7545        FORMAT("****** WARNING ! WARNING ! WARNING: " //
     &           " starting value for driven", i4,"-th bending is ",
     &           f10.4, " while actual value is ", f10.4 //
     &           "****** WARNING ! WARNING ! WARNING: " //)    
          END IF
        end do
      end if 
      if(nitors_added.gt.0) THEN 
        do i=1,nitors_added
          aux = abs(force_died(i)*(eqdied(i)-tors_added(i)*180./pi))
          if(aux.gt.400.d0) THEN 
            write(kprint,7546) i,eqdied(i),tors_added(i)*180./pi 
7546        FORMAT(//"****** WARNING ! WARNING ! WARNING: " // 
     &           " starting value for driven", i4,"-th torsion is ",
     &           f10.4, " while actual value is ", f10.4 // 
     &           "****** WARNING ! WARNING ! WARNING: " //)    
          END IF
        end do
      end if 

c---- check bonds

      if(steer_temperature) THEN 
        nsteer_temp=0
        write(kprint,1673) tiniz,tfina,timetiniz
     &       ,timetfina
        velt=(tfina-tiniz)/(timetfina-timetiniz)
        worktemp=0.d0 
        if(timetiniz.gt.timetfina) THEN 
          errmsg='Final time for non eq. thermal changes is smaller'
     &        //' than init time'
          CALL xerror(errmsg,80,1,2)
        end if
        if(t.ne.tiniz) THEN 
          errmsg='Current temperature not equal to inizial steer temp;'
     &        // ' program continues'
          CALL xerror(errmsg,80,1,21)
        end if
      end if

      if(kplot_stn0.gt.0) write(kplot_stn0,1011) 

c---  do input for program EE and DYNAMIC
 
      if(dynamic.and.slt_exist) THEN 
        if(lconstr.ne.0) THEN 
         errmsg=' Dynamic input not produced; constraints were found.'
         CALL xerror(errmsg,80,1,1)
       ELSE
         call do_dynamic_input (xp0,yp0,zp0,nato_slt,nmol_slv
     &        ,unitc,unite,avogad,kdynamic)
       END IF
      end if


*=======================================================================
*----- Choose on which shell to carry out the gmgp correction ----------
*=======================================================================

      h_skin=.FALSE.
      l_skin=.FALSE.
      m_skin=.TRUE.

*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot,protl)

*=======================================================================
*----- Change the cell according to input commands ---------------------
*=======================================================================

      IF(change_cell) THEN
         CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &        ,zpcma)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)

         CALL dcopy(9,co2,1,co,1)
         CALL dcopy(9,oc2,1,oc,1)

         CALL change_origin(-1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_frame(co,oc,1,nprot,xpcma,ypcma,zpcma,xpcm,ypcm
     &        ,zpcm)
      END IF

*=======================================================================
*---- Calculate group position  ----------------------------------------
*=======================================================================

      CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
      if(.not.scalemass)  CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass
     &     ,tmass,vcax,vcay,vcaz)
      
*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================

      CALL fndgrp(nprot,protl,atomp)

      IF(nflag(1) .NE. 0 .AND. ngrp_old .NE. ngrp ) THEN
         write(kprint,98000) ngrp_old,ngrp
      END IF
      IF(nflag(1) .NE. 0 .AND. nprot_old .NE. nprot ) THEN
         write(kprint,99000) nprot_old,nprot
      END IF

*=======================================================================
*----- If the system is at constant temperature initialize some --------
*----- variables and convert masses of the 3 thermostats ---------------
*=======================================================================

      IF(thermos) THEN
         CALL find_thermos(ntap,lconstr,natom_slt,natom_slv,nmol_slt
     &        ,nmol_slv,cnstpp_slv,nprot,cpress,isostress,ndf_thermos)
     &        
         CALL cov_thermos(slv_exist,slt_exist,qmass,ndf_thermos,t)
         IF(flag .EQ. 0) THEN
            CALL zero0(eta,neta)
         END IF
      END IF


*=======================================================================
*-------- Initialize velocities when starting without a restart --------
*=======================================================================

      IF(FLAG.EQ.0.AND.t.GT.1d-15) THEN
         LINIT=.TRUE.
         CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,linit
     &        ,massinfty)
         IF(cnstpp .NE. 0) THEN
            CALL rattle_correc(time,xp0,yp0,zp0,vpx,vpy
     &           ,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &           ,dnit,cnst_protp,cnst_protl,mim_lim,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF
         CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &        ,vcaz)
         LINIT=.FALSE.

*=======================================================================
*----- Set to zero velocities of the barostat --------------------------
*=======================================================================

         IF(cpress) THEN
            CALL set_tempp(masspr,vco,temppra,0.0D0)
         END IF
         IF(thermos) THEN
            CALL set_tempt(neta,qmass,etap,temph,0.0D0)
         END IF
      END IF   

*=======================================================================
*--- Compute the volume of the system ----------------------------------
*=======================================================================
      
      CALL matinv(3,3,co,oc,volume)
      volume=volume*boxl**3
     
*========================================================================
*==== Calls init routine for conventional kspace Ewald or PME -----------
*========================================================================

      IF(clewld) THEN

         IF(pme) THEN
            numatoms=ntap
#ifdef FFTW
            write(kprint,1003) 
            CALL fft_pme_init(numatoms,nfft1,nfft2,nfft3,pme_order
     &           ,bsp_mod1,bsp_mod2,bsp_mod3,planf,planr)
#else
            write(kprint,1004) 
            CALL fft_pme_init(numatoms,nfft1,nfft2,nfft3,pme_order
     &           ,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
     &           ,bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork)
#endif
            if ( siz_Q .GT. MAXT ) THEN
               write(kprint,78410)
               stop
            END IF
            rshk=shell_pme
            IF(erf_corr) THEN
              CALL erf_corr_cutoff(oc,delew,rkcut,nfft1,nfft2,nfft3)
              CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &             ,erf_arr_corr,work)
            END IF
         ELSE
            aux=(rtolh+rcuth)**2
            CALL chkewl(oc,aux,alphal,rkcut,volume)
            rshk=shell_pme
            IF(erf_corr) THEN  
              CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &           ,erf_arr_corr,work)
            END IF
         END IF

*=======================================================================
*----- Set up table for erfc spline ------------------------------------
*=======================================================================

         IF(erfc_spline) THEN
            aux=rcuth+rtolh
            aux=aux*aux
            CALL erfc_spline_init(aux,alphal,erfc_bin,mspline,erfc_arr
     &           ,work,iret,errmsg,rkcut,rcut_corr,erfc_spline_corr)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF
      END IF

*=======================================================================
*---- Calculate ewald self and the intramolecular terms. ---------------
*---- Only the contribution to the latter coming from  -----------------
*---- bond constraints is computed at this point -----------------------
*=======================================================================

      
      IF(clewld) THEN
         naus= 0
         if(alchemy) THEN
           do i=1,ntap
             chrge0(i)= chrge(i)*(1.d0 - dabs(lambdaq(i)))
           end do
         ELSE
           do i=1,ntap
             chrge0(i)= chrge(i)
           end do
         END IF
         CALL cself(ss_index,ntap,alphal,rkcut,chrge0,self_slt
     &        ,self_slv)
         CALL ferrf(ss_index,alphal,chrge0,1.0D0,xp0,yp0,zp0,0,lcnstr
     &        ,lconstr,fscnstr_slt,fscnstr_slv,fpx_p,fpy_p,fpz_p
     &        ,erf_corr,erf_arr_corr,delew,rlew,naus,ntap)
      END IF


*=======================================================================
*---  set fake neighbor lists
*=======================================================================

#ifdef _OMP_
      do i=1,nthr
        nnlppf(1,i)=ngrp
        DO j=1,ngrp
          nnlppf(1+j,i)=j
        end do
        do j=ngrp+2,mpp8
          nnlppf(j,i)=0
        end do
      end do
#else      
      nnlppf(1)=ngrp
      DO j=1,ngrp
         nnlppf(1+j)=j
      end do
#endif

*=======================================================================
*---  Check if there are no error in the list of bonds, bendings -------
*---  and torsions -----------------------------------------------------
*=======================================================================

      IF(prttopl) THEN
         CALL check_topology('B',top_bonds,lbond,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('D',top_bendings,lbend,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('P',top_ptors,ltors,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('I',top_itors,litor,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
      END IF


*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================

#ifdef _MPI_
! REM environment
      if(rem_run) then

        if( cdist.eq.0.OR.cdist.eq.3 ) then
          line_res = 'Simulation restarted from a previous one'
        else if  (cdist.eq.1.OR.cdist.eq.2) THEN  
          line_res = 'Simulation started from scratch'
        endif  
        ntraj_eff=ntrajectories/nbatteries
        ibatt=int((para_index-1)/ntraj_eff) + 1
        para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
        write(kprint,8000) iproc, nbatteries
     &       ,ntrajectories,ntraj_eff, t,t/rem_factor_max(1),t
     &       /rem_factor_max(2), t/rem_factor_max(3), rem_ts,rem_print
     &       ,rem_adjust_ratio,para_index,para0,ibatt,t/rem_factor(1),t
     &       /rem_factor(2), t/rem_factor(3), lskip_force, line_res
        if(cdist.eq.2) THEN 
          write(kprint,'(A57)')
     &       "-----  WARNING!! scale factors initializated for REM.set"
          write(kprint,'(A17,A57)')
     &         "-----            ", 
     &         "Active scale factors are those given in the REM.set!!!!"

        end if
        if( rem_groups >= 1 ) then
          write(kprint,8001)
          tlscaled=0
          do i = 1, rem_groups
            write(kprint,8002) i,rem_group(1,i),prsymb(nres(rem_group(1
     &           ,i),2)),rem_group(2,i),prsymb(nres(rem_group(2,i)
     &           ,2))
            tlscaled=tlscaled+rem_group(2,i)-rem_group(1,i)+1
          enddo
          write(kprint,'(8x,"scaled improper (and h-proper) torsions")') 
          do i=1,litor
            if(remitor(i)) THEN
              write(kprint,'(7x," - ",i5,4i5)') i,(litr(j,i),j=1,4)  
            endif
          end do
          write(kprint,'(/,7x," scaled (no-h) proper torsions")') 
          do i=1,ltors
            if(remltor(i)) THEN
              write(kprint,'(7x," - ",i5,4i5)') i,(ltor(j,i),j=1,4)  
            endif
          end do
          write(kprint,8003) tlscaled
          if( rem_segkind.eq.0) write(kprint,8004) 
          if( rem_segkind.eq.1) write(kprint,8005) 
          if( rem_segkind.eq.2) write(kprint,8006)
        endif
      endif
#endif


! SGE environment
      if ( sim_tempering ) THEN 
        if(always_accept) THEN 
          write(kprint,8007) target_replica,2
     &         *ntrajectories*sge_ts
          nplot_fragm=0
          nplot_center=0
          nascii=0
        ELSE
          call sge_print_titles ( iproc, ntrajectories,
     &         atom_b1, atom_b2, atom_be1, atom_be2, atom_be3,
     &         atom_it1, atom_it2, atom_it3, atom_it4 )
        END IF
      end if
      
! Metadynamics environment
      IF(meta_run) THEN
        IF (meta_analysis) THEN
          line_res = 'Metadynamics: analysis run'
        ELSE
          IF(metar) THEN
            line_res =
     &           'Metadynamics simulation restarted from a previous one'
          ELSE
            line_res =
     &           'Metadynamics simulation started from scratch'
          END IF
        END IF
        write(kprint,9000) line_res, meta_nrc, meta_ts,
     &       meta_height *(unite*avogad) / 1000.0,(meta_width(j), j = 1
     &       ,meta_nrc)
      END IF
      



      CALL prtite
      IF(prttpg) THEN
         write(kprint,1000)
         write(kprint,'(////)')
         IF(prtseq) CALL prtsq(nbun,mend,prsymb)
         IF(prtatl) CALL prtat(ss_index,ntap,beta,betb,nres,m1,prsymb
     &        ,chrge,mass)
         IF(prtbndl) CALL prtbnd(beta,nres(1,1),lstrtch,lstretch,potbo
     &        ,m9)
         IF(prtcnl) CALL prtcn(beta,nres(1,1),lcnstr,lconstr,potbo_cnst
     &        ,m9)
         IF(prtbal) CALL prtba(beta,nres(1,1),lbndg,lbend,potbe,m2)
         IF(prtptl) CALL prtpt(beta,nres(1,1),ltor,ltors,potto,m3)
         IF(prtitl) CALL prtit(beta,nres(1,1),litr,litor,potit,m4)
         write(kprint,'(////77(''=''))')
      END IF

*======================================================================
*==== Write force field when required. Only the parameters used -------
*==== in the actual simulation will be printed ------------------------
*======================================================================

      IF(write_ff_pars) THEN
         write(kprint,89000)
         write(kprint,*) "--> Begin PRM"
         write(kprint,90000)
         call write_pot_bond(lstretch,2,lstrtch,m9,2
     &        ,potbo,betb,1)
         write(kprint,95000)
         call write_pot_bond(lconstr,2,lcnstr,m9,2
     &        ,potbo_cnst,betb,-1)
         write(kprint,94000)
         write(kprint,91000)
         call write_pot_bond(lbend,3,lbndg,m2,4,potbe
     &        ,betb,1)
         write(kprint,94000)
         write(kprint,92000)
         call write_pot_bond(ltors,4,ltor,m3,2,potto
     &        ,betb,1)
         write(kprint,94000)
         write(kprint,93000)
         call write_pot_bond(litor,4,litr,m4,3,potit
     &        ,betb,1)
         write(kprint,94000)
         write(kprint,96000)
         write(kprint,97000)
         call write_pot_nbond(ecc6,ecc12,ecc146,ecc1412
     &        ,nbtype,mass,lj_fudgeb,betb,ntap)
         write(kprint,94000)
         write(kprint,*) "--> End PRM"
         write(kprint,89500)
      END IF

*======================================================================
*==== Write a Banner Before intermediate results ----------------------
*======================================================================

      write(kprint,'(//)')
      write(kprint,1200)
      write(kprint,1300)
      write(kprint,1110)
      write(kprint,1110)
      write(kprint,1300)
      write(kprint,1200)
      write(kprint,'(//)')


      IF(rneih.GT.0.d0) THEN
         lupdate=.true.
      else
         nnlpp0(1,1)=ngrp
         DO j=1,ngrp
            nnlpp0(1+j,1)=j
         end do
         lupdate=.false.
      END IF   

*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)

*=======================================================================
*--- Compute lx, ly, lz by subtracting the coordinates of the c.of m.---
*=======================================================================

      CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma,ypcma
     &     ,zpcma,co)

*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE)          

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse
      IF(lupdate.and.(.not.linked_cell)) THEN
*--      update shell h neighbor list
        CALL mts_forces(nstep,'u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &       ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &       ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p,stressd_p
     &       ,nnlppf,nnlpp0)
      ELSE
         aux= rcuth+rtolh+rneih
#ifdef _OMP_
         npoints=0
!        call cell mapping 
         call mapp_lc(ncx,ncy,ncz,i_ind,j_ind,k_ind,n_ind,npoints)
         CALL lc_index_omp(indmax,ncx,ncy,ncz,nind,indxi,indxj,indxk
     &        ,i_ind,j_ind,k_ind,n_ind,aux,co,npoints)
#else
         CALL lc_index(indmax,ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co)
#endif
         CALL timer(vfcp,tfcp,elapse)
         gcpu=-gcpu + tfcp
         gela=-gela + elapse
         write(kprint,15011) gcpu,gela
         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         gela=elapse
         CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &        ,ypga,zpga,ngrp,nnlpp0,kprint)
      END IF   

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gela=-gela + elapse
      gcpu_u=gcpu
      gela_u=gela
      write(kprint,16011) gcpu,gela_u

*=======================================================================
*---- Zeroes all forces ------------------------------------------------
*=======================================================================

      CALL zeroa(fpx_h,fpy_h,fpz_h,ntap,1)
      CALL zeroa(fpx_l,fpy_l,fpz_l,ntap,1)
      CALL zeroa(fpx_m,fpy_m,fpz_m,ntap,1)
      CALL zeroa(fpx_n1,fpy_n1,fpz_n1,ntap,1)
      CALL zeroa(fpx_n0,fpy_n0,fpz_n0,ntap,1)

      write(kprint,15000) 

*=======================================================================
*---- Do check on coordinates before starting  if required
*=======================================================================

      if(check_coord) THEN
        call check_coordinates(ntap,co,xpa,ypa,zpa,imin,jmin
     &       ,laux,drmin)
        write(kprint,40098) 
40098   format(/79("-"),/,20(" "), "Coordinate check",/79("-"))    
        if(laux.and.(.not.alchemy)) THEN 
          write(kprint,40909) imin,jmin, dsqrt(drmin)  
40909    FORMAT( 10x," Atom", i8, " and", i8,
     &           " are too close: r = ",E12.5)
          errmsg='Problems with coordinates; Program may fail... '
          CALL xerror(errmsg,80,2,11)
          errmsg='!!! Atomic coordinates written on fort.888 file !!'
          CALL xerror(errmsg,80,2,11)
          do i=1,ntap
            j = nres(i,2)
            write(888,31003)'ATOM  ',i,beta(i)(1:5),prsymb(j)(1:3)
     &           ,nres(i,1),xp0(i),yp0(i),zp0(i)
31003       FORMAT(a6,i5,1x,a5,a3,2x,i4,4x,3f8.3)
          end do
        end if
        close(888)
        write(kprint,40099) 
40099   format(/79("-"))

        if(.not.laux) THEN 
5008      FORMAT(80A1) 
          write(kprint,6754) sqrt(drmin),imin,jmin 
6754      FORMAT(80("="),/"=",78x,"=",/,
     &         "=  Initial coordinates OK! ",52x,"=",/
     &         "=  Shortest distance   ",f10.3," between atoms",i6,1x,i6
     &         ,19x,"=",/,"=",78x,"=",/80("="))
        END IF
      ENDIF


*=======================================================================
*===== Computes all forces at time = 0 ================================= 
*=======================================================================

************************************************************************
***                        Direct SPACE                              ***
***        keep this order to compute nested neighbor lists          ***
************************************************************************

*=======================================================================
*==== H-contribution in direct space
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

*---  CALLs mts_forslv,mts_fnbond,mts_forpw

      
      rshell='h'
      CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_h,ucos_h,virs_h,virsp_h
     &     ,ucnp_h,ucop_h,ucnsp_h,ucosp_h,fpx_h,fpy_h,fpz_h,stressd_h
     &     ,nnlpp0,nnlpp)
      CALL timer(vfcp,tfcp,elapse)
      gcpu_hd=-gcpu + tfcp
      gcpu_1 = gcpu_hd
      gela_hd=-gela + elapse
      gela_1 = gela_hd

*=======================================================================
*==== L-contribution in direct space (long-ranged) ---------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

      rshell='l'
      CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_l,ucos_l,virs_l,virsp_l
     &     ,ucnp_l,ucop_l,ucnsp_l,ucosp_l,fpx_l,fpy_l,fpz_l,stressd_l
     &     ,nnlpp,nnlpp1)
      CALL timer(vfcp,tfcp,elapse)
      gcpu_ld=-gcpu + tfcp
      gcpu_2 = gcpu_ld
      gela_ld=-gela + elapse
      gela_2 = gela_ld

*=======================================================================
*------- M-contribution in direct space (short-ranged) -----------------
*-------          and protein torsion ----------------------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse
      rshell='m'
      CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_m,ucos_m,virs_m,virsp_m
     &     ,ucnp_m,ucop_m,ucnsp_m,ucosp_m,fpx_m,fpy_m,fpz_m,stressd_m
     &     ,nnlpp1,nnlpp2)

      IF(abmd .AND. (dissociate .OR. associate)) CALL
     &     comp_abmd_fdiss(abmd_dir,rspset,diss_list,spring,xpa,ypa,zpa
     &     ,co,fpx_m,fpy_m,fpz_m,uumb,gr)

      CALL timer(vfcp,tfcp,elapse)
      gcpu_md=-gcpu + tfcp
      gcpu_3 = gcpu_md
      gela_md=-gela + elapse
      gela_3 = gela_md

************************************************************************
***                    Reciprocal SPACE                              ***
***     order is reverted to avoid to compute twice pme corrections  ***
************************************************************************

*=======================================================================
*==== M-contribution in reciprocal Space -------------------------------
*=======================================================================

      IF(.NOT.clewld) goto 907 
      gcpu=0.d0
      gela=0.d0
      rshell='m'
      IF(pme .AND. rshell .NE. rshk) goto 914
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

*---  CALLs mts_furier, mts_furipp, mts_furipw 


      if(alchemy) THEN
        CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &       ,zpcma,urcsp_m,urcs_m,urcp_m,virsp_m,virs_m,virp_m,fpx_m
     &       ,fpy_m,fpz_m,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt_rl,coul_bnd_slv
     &       ,rshell,rshk,eer_m,stressr_m,fudgec,self_slt)
      ELSE
        CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &       ,zpcma,urcsp_m,urcs_m,urcp_m,virsp_m,virs_m,virp_m,fpx_m
     &       ,fpy_m,fpz_m,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell
     &       ,rshk,eer_m,stressr_m,fudgec)
      END IF
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gela=-gela + elapse
914   write(kprint,16003) gcpu,gcpu_md,gcpu+gcpu_md,gela,gela_md,gela
     &     +gela_md
      gcpu_1=gcpu+gcpu_md
      gela_1=gela+gela_md

*=======================================================================
*==== L-contribution in reciprocal Space -------------------------------
*=======================================================================

      gcpu=0.d0
      gela=0.d0
      rshell='l'
      IF(pme .AND. rshell .NE. rshk) goto 915
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

*--   subtract the less accurate k-forces

      IF(alchemy) THEN 
        CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &       ,zpcma,urcsp_l,urcs_l,urcp_l,virsp_l,virs_l,virp_l,fpx_l
     &       ,fpy_l,fpz_l,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt_rl,coul_bnd_slv
     &       ,rshell,rshk,eer_l,stressr_l,fudgec,self_slt)
      ELSE
        CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &       ,zpcma,urcsp_l,urcs_l,urcp_l,virsp_l,virs_l,virp_l,fpx_l
     &       ,fpy_l,fpz_l,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell
     &       ,rshk,eer_l,stressr_l,fudgec)
      END IF
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gela=-gela + elapse
915   write(kprint,16002) gcpu,gcpu_ld,gcpu+gcpu_ld,gela,gela_ld,gela
     &     +gela_ld
      gcpu_2=gcpu+gcpu_ld
      gela_2=gela+gela_ld


*=======================================================================
*==== H-contribution in reciprocal Space -------------------------------
*=======================================================================

      gcpu=0.d0      
      gela=0.d0
      rshell='h'
      IF(pme .AND. (rshell .NE. rshk)) goto 916
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

*--   subtract the less accurate k-forces

      if(alchemy) THEN 
        CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &       ,zpcma,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h
     &       ,fpy_h,fpz_h,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt_rl,coul_bnd_slv
     &       ,rshell,rshk,eer_h,stressr_h,fudgec,self_slt)
      ELSE
        CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma
     &       ,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h,fpy_h
     &       ,fpz_h,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &       ,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell
     &       ,rshk,eer_h,stressr_h,fudgec)
      END IF

c---  if alchemy computes starting work  
      CALL timer(vfcp,tfcp,elapse)

      gcpu=-gcpu + tfcp
      gela=-gela + elapse
916   write(kprint,16001) gcpu,gcpu_hd,gcpu+gcpu_hd,gela,gela_hd,gela
     &     +gela_hd
      gcpu_3=gcpu+gcpu_hd
      gela_3=gela+gela_hd
907   theoric_speed_up=(gela_1+gela_2+gela_3)*mrespa*lrespa/(gela_1
     &     *mrespa*lrespa+gela_2*lrespa+gela_3)

#ifdef _MPI_
      if( rem_run .or. (sim_tempering .and. .not.sim_tempering_sge) ) 
     &     then
#else
      if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
        if(.not.rem_segment) then 
!if rem_segment is .TRUE. the real part of the nb forces are scaled
!inside mts_forpp_rem, while the reciprocal contribution is not scaled
          call dscal(ntap,rem_factor(3),fpx_h,1)
          call dscal(ntap,rem_factor(3),fpy_h,1)
          call dscal(ntap,rem_factor(3),fpz_h,1)
          call dscal(ntap,rem_factor(3),fpx_l,1)
          call dscal(ntap,rem_factor(3),fpy_l,1)
          call dscal(ntap,rem_factor(3),fpz_l,1)
          call dscal(ntap,rem_factor(3),fpx_m,1)
          call dscal(ntap,rem_factor(3),fpy_m,1)
          call dscal(ntap,rem_factor(3),fpz_m,1)
          do i = 1,3
            do j = 1,3
              stressd_h(i,j) = stressd_h(i,j)*rem_factor(3)
              stressr_h(i,j) = stressr_h(i,j)*rem_factor(3)
              stressd_l(i,j) = stressd_l(i,j)*rem_factor(3)
              stressr_l(i,j) = stressr_l(i,j)*rem_factor(3)
              stressd_m(i,j) = stressd_m(i,j)*rem_factor(3)
              stressr_m(i,j) = stressr_m(i,j)*rem_factor(3)
            enddo
          enddo
        endif
      endif
      
*=======================================================================
*--- Compute force on the co matrix ------------------------------------
*=======================================================================

      IF(cpress) THEN
         CALL comp_stress_conf(stressd_m,stressr_m,prt_m,oc,volume,unitp
     &        ,press_m)
         CALL comp_stress_kinetic(vcax,vcay,vcaz,tmass,co,nprot,volume
     &        ,unitp,st_m,press_kin)
         CALL comp_forcep(prt_m,st_m,oc,volume,pext)
         CALL comp_stress_conf(stressd_l,stressr_l,prt_l,oc,volume,unitp
     &        ,press_l)
         CALL comp_stress_conf(stressd_h,stressr_h,prt_h,oc,volume,unitp
     &        ,press_h)
      END IF

*=======================================================================
*--- Print timing for spherical cutoff ---------------------------------
*=======================================================================

      IF(.NOT.clewld) THEN
         gcpu=0.d0
         gela=0.d0
         write(kprint,16001) gcpu,gcpu_hd,gcpu_hd,gela,gela_hd,gela_hd
         write(kprint,16002) gcpu,gcpu_ld,gcpu_ld,gela,gela_ld,gela_ld
         write(kprint,16003) gcpu,gcpu_md,gcpu_md,gela,gela_md,gela_md
         gcpu_1=gcpu+gcpu_md
         gcpu_2=gcpu+gcpu_ld
         gcpu_3=gcpu+gcpu_hd
         gela_1=gela+gela_md
         gela_2=gela+gela_ld
         gela_3=gela+gela_hd
      END IF   
      write(kprint,10067) theoric_speed_up

************************************************************************
***            FAST INTRAMOLECULAR COMPONENTS                        ***
************************************************************************

*=======================================================================
*---- Computes bonded forces of shell n1 -------------------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

      CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n1,fpy_n1
     &     ,fpz_n1,fudge,lj_fudge,abmd_dir,puhyd,conf_bnd_slt_n1
     &     ,conf_bnd_slv_n1,coul_bnd_slt_n1,coul_bnd_slv_n1,unb14,cnb14
     &     ,ungrp,cngrp,uptors,uslvtor,mapdn,nmapdn,uumb,gr)
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gela=-gela + elapse
      write(kprint,14004) gcpu,gela
      gcpu_0=gcpu
      gela_0=gela

#ifdef _MPI_
      if( rem_run .or. (sim_tempering .and. .not.sim_tempering_sge) ) 
     &     then
#else
      if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
        if(.not.rem_segment) then  !if rem_segment is .TRUE. the n1 forces are scaled inside fnb14_rem and fptor_rem
          call dscal(ntap,rem_factor(2),fpx_n1,1)
          call dscal(ntap,rem_factor(2),fpy_n1,1)
          call dscal(ntap,rem_factor(2),fpz_n1,1)
        endif
      endif

*=======================================================================
*---- Computes bonded forces of shell n0 -------------------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      gela=elapse

      CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n0,fpy_n0
     &     ,fpz_n0,ubond,uslvbon,ubend,uslvben,uitors,uslvitor,ndim)

      if(addstrcom) THEN 
        laux=.true.  
        call com_restraint (kcom_md,r0com,xp0,yp0,zp0,mlig,mtar
     &       ,mass,rcom,fpx_n0,fpy_n0,fpz_n0,ucom,ilig1
     &       ,ilig2,itar1,itar2,laux)
      ENDIF

      CALL timer(vfcp,tfcp,elapse)

#ifdef _MPI_
      if( rem_run .or. (sim_tempering .and. .not.sim_tempering_sge) ) 
     &     then
#else
      if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
        if  (.not.rem_segment) THEN 
          call dscal(ntap,rem_factor(1),fpx_n0,1)
          call dscal(ntap,rem_factor(1),fpy_n0,1)
          call dscal(ntap,rem_factor(1),fpz_n0,1)
        END IF
      endif

*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for all shell forces --------------------------
*=======================================================================

      CALL comp_fcm(nprot,protl,fpx_m,fpy_m,fpz_m,fcax_m,fcay_m,fcaz_m
     &     ,mass,tmass,oc)

      CALL comp_fcm(nprot,protl,fpx_l,fpy_l,fpz_l,fcax_l,fcay_l,fcaz_l
     &     ,mass,tmass,oc)

      CALL comp_fcm(nprot,protl,fpx_h,fpy_h,fpz_h,fcax_h,fcay_h,fcaz_h
     &     ,mass,tmass,oc)

!     this is done since added topology may couple different molecule
!     so that there is a net force on the COM due to the added topol
      IF(nadd_tpg.gt.0.or.addstrcom) THEN
        CALL comp_fcm(nprot,protl,fpx_n0,fpy_n0,fpz_n0,fcax_n0,fcay_n0
     &        ,fcaz_n0,mass,tmass,oc)
      END IF

*=======================================================================
*---- Print timing and simulation CPU time estimates -------------------
*=======================================================================
      
      gcpu=-gcpu + tfcp
      gela=-gela + elapse
      write(kprint,16004) gcpu,gela
      theoric_speed_up=(gela+gela_0+gela_1+gela_2+gela_3)*n0respa
     &     *n1respa*mrespa*lrespa/(gela*n0respa*n1respa*mrespa*lrespa
     &     +gela_0*n1respa*mrespa*lrespa+gela_1*mrespa*lrespa+gela_2
     &     *lrespa+gela_3)
      write(kprint,10068) theoric_speed_up
      if(nupdte.gt.0) gela_u = gela_u /dfloat(nupdte)
      aux = (gela*n0respa*n1respa*lrespa*mrespa+gela_0*n1respa*lrespa
     &     *mrespa+gela_1*mrespa*lrespa+gela_2*lrespa + gela_3 + gela_u)
      gela_0 = aux/dfloat(lrespa*mrespa)
      gela_1 = aux / time
      aux = aux*(maxstp+(mrject/dfloat(lrespa*mrespa)))
      hours = int(aux/3600.d0)
      min = int((aux-hours*3600)/60.d0)
      write(kprint,80033) hours,min,gela_0,gela_1
      write(kprint,60030) 

      CALL timer(vfcp,tfcp,elapse)

      if(sim_tempering.AND. always_accept) THEN 
         CALL dcopy(9,co,1,co_s,1)
         CALL dcopy(9,oc,1,oc_s,1)
         CALL dcopy(ntap,xp0,1,xp0_s,1)
         CALL dcopy(ntap,yp0,1,yp0_s,1)
         CALL dcopy(ntap,zp0,1,zp0_s,1)
         CALL dcopy(ntap,lx,1,lx_s,1)
         CALL dcopy(ntap,ly,1,ly_s,1)
         CALL dcopy(ntap,lz,1,lz_s,1)
         CALL dcopy(ngrp,xpg,1,xpg_s,1)
         CALL dcopy(ngrp,ypg,1,ypg_s,1)
         CALL dcopy(ngrp,zpg,1,zpg_s,1)
         CALL dcopy(nprot,xpcma,1,xpcma_s,1)
         CALL dcopy(nprot,ypcma,1,ypcma_s,1)
         CALL dcopy(nprot,zpcma,1,zpcma_s,1)
         CALL dcopy(ntap,fpx_m,1,fpx_m_s,1)
         CALL dcopy(ntap,fpy_m,1,fpy_m_s,1)
         CALL dcopy(ntap,fpz_m,1,fpz_m_s,1)
         CALL dcopy(ntap,fpx_l,1,fpx_l_s,1)
         CALL dcopy(ntap,fpy_l,1,fpy_l_s,1)
         CALL dcopy(ntap,fpz_l,1,fpz_l_s,1)
         CALL dcopy(ntap,fpx_h,1,fpx_h_s,1)
         CALL dcopy(ntap,fpy_h,1,fpy_h_s,1)
         CALL dcopy(ntap,fpz_h,1,fpz_h_s,1)
         CALL dcopy(ntap,fpx_n1,1,fpx_n1_s,1)
         CALL dcopy(ntap,fpy_n1,1,fpy_n1_s,1)
         CALL dcopy(ntap,fpz_n1,1,fpz_n1_s,1)
         CALL dcopy(ntap,fpx_n0,1,fpx_n0_s,1)
         CALL dcopy(ntap,fpy_n0,1,fpy_n0_s,1)
         CALL dcopy(ntap,fpz_n0,1,fpz_n0_s,1)
         CALL dcopy(9,stressd,1,stressd_s,1)
         CALL dcopy(9,stressr,1,stressr_s,1)
         CALL dcopy(9,gmgp,1,gmgp_s,1)
         CALL dcopy(3,eta,1,eta_s,1)
         CALL dcopy(3,fth,1,fth_s,1)
         CALL dcopy(1,temp,1,temp_s,1)
       END IF
*=======================================================================
*===== MD loop starts here =============================================
*=======================================================================
      gcpu=tfcp
      elaps=elapse
      lfirst=.true.
      time2=time*0.5D0

      if(energy_then_die) THEN 
        time=0.d0 
        nprint=1
      end if
        
100   CONTINUE


      nstep=nstep+1
      mstep=mstep+1

*=======================================================================
*---- Start 4 time steps r-RESPA integrator ----------------------------
*=======================================================================

*=======================================================================
*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE) --
*=======================================================================

      IF(MOD(nstep,nupdte) .EQ. 0 .AND. lupdate) THEN
         
*=======================================================================
*--  1)  Update shell *H* neighbor list
*=======================================================================
         
         IF(.not.linked_cell) THEN
           CALL mts_forces(nstep,'u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &          ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &          ,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &          ,stressd_p,nnlppf,nnlpp0)
            rcuth_save=rcuth
            rtolh_save=rtolh
            rcuth=rcutl
            rtolh=0.0D0
         ELSE   
            aux=rcuth+rtolh+rneih 
            if(omp_timing) THEN 
              call timer(treal,tcpu1,tela1)
            end if
            IF(cpress .AND. MOD(nstep,nupdte*nupdte_index) .EQ. 0) THEN
#ifdef _OMP_
              CALL lc_index_omp(indmax,ncx,ncy,ncz,nind,indxi,indxj
     &             ,indxk,i_ind,j_ind,k_ind,n_ind,aux,co,npoints)
              
#else
               CALL lc_index(indmax,ncx,ncy,ncz,nind,indxi,indxj,indxk
     &              ,aux,co)
#endif
               write(kprint,2000)
            END IF
            CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &           ,ypga,zpga,ngrp,nnlpp0,kprint)
            if(omp_timing) THEN 
              call timer(treal,tcpu2,tela2)
              time_u=time_u+tela2-tela1
            end if
            rcuth_save=rcuth
            rtolh_save=rtolh
            rcuth=rcutl
            rtolh=0.0D0
         END IF 
         
*=======================================================================
*--  2)  Update shell *L* neighbor list
*=======================================================================
         
         CALL mts_forces(nstep,'h',xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &        ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &        ,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &        ,stressd_p,nnlpp0,nnlpp)
         rcutl_save=rcutl
         rtoll_save=rtoll
         rcutl=rcutm
         rtoll=0.0D0
         
*=======================================================================
*--  3)  Update shell *M* neighbor list
*=======================================================================
         
         CALL mts_forces(nstep,'l',xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &        ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &        ,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &        ,stressd_p,nnlpp,nnlpp1)
         rcuth=rcuth_save
         rtolh=rtolh_save
         rcutl=rcutl_save
         rtoll=rtoll_save
      END IF   
      
*=======================================================================
*---  Advances velocities for half time step TIME using H-forces    ----
*=======================================================================
      
      CALL correc(vpx,vpy,vpz,fpx_h,fpy_h,fpz_h,mass,ntap,time)
      
      IF(cnstpp .NE. 0) THEN
         CALL rattle_correc(time,xp0,yp0,zp0,vpx,vpy
     &        ,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass,dnit
     &        ,cnst_protp,cnst_protl,mim_lim,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF

      CALL correc_stress(cpress,h_skin,nprot,co,oc,vco,vcax,vcay
     &     ,vcaz,fcax_h,fcay_h,fcaz_h,stressd_h,stressr_h,volume
     &     ,press_h,press_kin,pext,tmass,masspp,time,time2)

      IF(isostress) THEN
        CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,nboxcnst,iret
     &       ,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF
      
*=======================================================================
*---  Reset h-forces to ZERO                                        ----
*=======================================================================
      
#ifdef _OMP_
      CALL zeroomp(fpx_h,fpy_h,fpz_h,ntap)
#else
      CALL zeroa(fpx_h,fpy_h,fpz_h,ntap,1)
#endif
      tl = time/dfloat(lrespa)
      tl2=tl*0.5D0
      
*=======================================================================
*---     Advances velocities for half time step *time* using        ----
*--- |-->         long range L-forces                               ----
*=======================================================================
      
      DO il=1,lrespa
        CALL correc(vpx,vpy,vpz,fpx_l,fpy_l,fpz_l,mass,ntap,tl)
         IF(cnstpp .NE. 0) THEN
            CALL rattle_correc(tl,xp0,yp0,zp0,vpx,vpy,vpz,ntap,cnstp
     &           ,dssp,coeffp,cnstpp,mass,dnit,cnst_protp,cnst_protl
     &           ,mim_lim,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

         CALL correc_stress(cpress,l_skin,nprot,co,oc,vco,vcax,vcay
     &        ,vcaz,fcax_l,fcay_l,fcaz_l,stressd_l,stressr_l,volume
     &        ,press_l,press_kin,pext,tmass,masspp,tl,tl2)

         IF(isostress) THEN
            CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,nboxcnst
     &           ,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

*=======================================================================
*---     Reset l-forces to ZERO                                    -----
*=======================================================================
         
#ifdef _OMP_
         CALL zeroomp(fpx_l,fpy_l,fpz_l,ntap) 
#else
         CALL zeroa(fpx_l,fpy_l,fpz_l,ntap,1) 
#endif
         tm = tl/dfloat(mrespa)
         tm2=tm*0.5D0
         tm4=tm*0.25D0
         
         rshell='m'
         
*=======================================================================
*---        Advances velocities for half time step *time* using     ----
*---    |-->           short range M-forces                         ----
*=======================================================================
         
         DO im=1,mrespa
            ninner=ninner + 1
            nsstt = nsstt + 1
            
            IF(thermos) THEN
               IF(isostress) THEN
                 CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &                ,nboxcnst,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap
     &                 ,cnstp,dssp,coeffp,cnstpp,mass,dnit,cnst_protp
     &                 ,cnst_protl,mim_lim,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
c---------- Correct velocities of the NETA thermostats ----------------
c---------- Propagate exp(i L_y tm/4) ---------------------------------
               CALL comp_thermos_forces(cpress,ndf_thermos,ntap,ss_index
     &              ,nprot,co,vpx,vpy,vpz,vcax,vcay,vcaz,vco,mass,tmass
     &              ,masspr,t,fth)
               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF
*----------  Propagate exp(i L_z tm/2) ---------------------------------
            IF(cpress .AND. coupl_mol) THEN
               CALL correc_matr(tm2,co,oc,vco,gmgp,vcax,vcay,vcaz,nprot)
            END IF
            IF(thermos) THEN
c----------  Correct velocities of the NETA thermostats ----------------
c----------  Recompute forces on the thermostat variables --------------
               CALL comp_thermos_forces(cpress,ndf_thermos,ntap,ss_index
     &              ,nprot,co,vpx,vpy,vpz,vcax,vcay,vcaz,vco,mass,tmass
     &              ,masspr,t,fth)
*==========  Propagate exp(i L_y tm/4) =================================
               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF
            
            IF(thermos) THEN
*=========  Propagate exp(i L_x tm/4) ==================================
               CALL correc_exp_scale(cpress,slt_exist,slv_exist
     &              ,ss_point(1,1),ss_point(1,2),etap,tm4,nprot,vcax
     &              ,vcay,vcaz,vpx,vpy,vpz,vco)
            END IF

*========= Propagate exp(i L_s + L_u tm/2) -==================================

            CALL correc(vpx,vpy,vpz,fpx_m,fpy_m,fpz_m,mass,ntap,tm)
            
            IF(cnstpp .NE. 0) THEN
              CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap,cnstp
     &              ,dssp,coeffp,cnstpp,mass,dnit,cnst_protp,cnst_protl
     &              ,mim_lim,iret,errmsg)
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF

            CALL correc_stress(cpress,m_skin,nprot,co,oc,vco,vcax,vcay
     &           ,vcaz,fcax_m,fcay_m,fcaz_m,stressd_m,stressr_m,volume
     &           ,press_m,press_kin,pext,tmass,masspp,tm,tm2)
            IF(isostress) THEN
              CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,nboxcnst
     &             ,iret,errmsg)
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF
            
*=========  Propagate exp(i L_x tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(cpress,slt_exist,slv_exist
     &              ,ss_point(1,1),ss_point(1,2),etap,tm4,nprot,vcax
     &              ,vcay,vcaz,vpx,vpy,vpz,vco)
            END IF
            
*=======================================================================
*---        Reset m-forces to ZERO
*=======================================================================
            
#ifdef _OMP_
            CALL zeroomp(fpx_m,fpy_m,fpz_m,ntap)
#else
            CALL zeroa(fpx_m,fpy_m,fpz_m,ntap,1)
#endif
            tn1 = tm/dfloat(n1respa)
            tn12=tn1*0.5D0
            
*=======================================================================
*---                Advance position and velocities for a full      ----
*---       |-->     time step tn1 using intramolecular forces       ----
*=======================================================================
            
            if(alchemy) dwrk14t=0.d0
            DO in1=1,n1respa
              ninn1=ninn1+1
              CALL correc(vpx,vpy,vpz,fpx_n1,fpy_n1,fpz_n1,mass,ntap
     &             ,tn1)
              IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(tn1,lx,ly,lz,vpx,vpy,vpz
     &                 ,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                 ,dnit,cnst_protp,cnst_protl,mim_lim,iret
     &                 ,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
#ifdef _OMP_
               CALL zeroomp(fpx_n1,fpy_n1,fpz_n1,ntap)
#else
               CALL zeroa(fpx_n1,fpy_n1,fpz_n1,ntap,1)
#endif
               tn0 = tn1/dfloat(n0respa)
               tn02=tn0*0.5D0
               
               DO in0=1,n0respa
                  ninn0=ninn0+1

#ifdef _OMP_
                  call copyomp(ntap,lx,ly,lz,xpo,ypo,zpo)
#else
                  CALL dcopy(ntap,lx,1,xpo,1) 
                  CALL dcopy(ntap,ly,1,ypo,1) 
                  CALL dcopy(ntap,lz,1,zpo,1) 
#endif                  
                  if(.not.start_conf)  then
                    CALL verlet(mass,ntap,lx,ly,lz,vpx,vpy,vpz
     &                   ,fpx_n0,fpy_n0,fpz_n0,tn0)
                  else
                    CALL starting_verlet(mass,ntap,lx,ly,lz,vpx,vpy
     &                    ,vpz,fpx_n0,fpy_n0,fpz_n0,tn0,max_dist)
                  endif
                  IF(cnstpp .NE. 0) THEN
                     CALL rattle_verlet(tn0,lx,ly,lz,xpo,ypo,zpo,vpx
     &                    ,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                    ,dnit,cnst_protp,cnst_protl,mim_lim,iret
     &                    ,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF

                  if(nadd_tpg.gt.0.or.addstrcom) THEN
                    CALL correc(vcax,vcay,vcaz,fcax_n0,fcay_n0,fcaz_n0
     &                   ,tmass,nprot,tn0)
                    if(cpress) call comp_stress_add_tpg(nprot,stresst0
     &                   ,co,oc,xpcm,ypcm,zpcm,fcax_n0,fcay_n0,fcaz_n0
     &                   ,vco,masspp,tn0)
                  ENDIF 

                  IF(isostress) THEN
                    CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &                   ,nboxcnst,iret,errmsg)
                    IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF

*----------  Propagate exp(i L_y tm/2) ---------------------------------

                  CALL verlet_free(nprot,xpcma,ypcma,zpcma,vcax,vcay
     &                 ,vcaz,tn0)

*---              Reset n-forces to ZERO
                  
#ifdef _OMP_
                  CALL zeroomp(fpx_n0,fpy_n0,fpz_n0,ntap) 
#else
                  CALL zeroa(fpx_n0,fpy_n0,fpz_n0,ntap,1) 
#endif                  
*=======================================================================
*---           Compute N0-forces                                    ----
*=======================================================================
                  
                  IF(cpress) THEN
                     DO i=1,3
                        DO j=1,3
                           coo(i,j)=co(i,j)
                        END DO
                     END DO
                     
*--- Propagate co matrix for step tm -----------------------------------
                     
                     CALL verlet_free(3,co(1,1),co(1,2),co(1,3),vco(1,1)
     &                    ,vco(1,2),vco(1,3),tn0)
                     IF(isostress) THEN
                        CALL rattle_verlet_co(tm,co,coo,dssco,cnstco,vco
     &                       ,masspp,nboxcnst,iret,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                     
                     CALL matinv(3,3,co,oc,volume)
                     volume=volume*boxl**3
!                     write(6,566) in0,in1,im,il,il,volume
!566                  FORMAT(5I5,G15.5)
                  END IF

*=======================================================================
*---  Recompute coordinates before recomputing forces ------------------
*=======================================================================

#ifdef _OMP_
                  CALL coord_inner_omp(ntap,nprot,iprot,co
     &                 ,oc,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa,xpcm,ypcm
     &                 ,zpcm,xpcma,ypcma,zpcma)                  
#else
                  CALL change_coord_inner(ntap,nprot,protl,ngrp,grppt,co
     &                 ,oc,pmass,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa,xpcm
     &                 ,ypcm,zpcm,xpcma,ypcma,zpcma) 
#endif
                  rtime= oldtime + (ninn0*tn0)

*--------         steered MD implemented (Jarzynsky on extra bonds bends
*--------         tors) 
                  if(time_steer) THEN
                    aux=1.d0
c========           standard steering on a linear path in N-dim
c========           reactions n coordinates
                    if(.not.path_steer.or.alchemy) THEN 
                      call steer_dynamics(lstretch,nbonds_added,strbonds
     &                     ,potbo(1,2),eqdist,eqdist1,str_vel,rtime
     &                     ,str_time0,str_time1,str_tim0r,str_tim1r,aux)
*--------             steered MD implemented (Jarzynsky on extra bends) 
                      call steer_dynamics(lbend,nbends_added,strbends
     &                     ,potbe(1,2),eqang,eqang1,ben_vel,rtime
     &                     ,str_time0,str_time1,str_tim0r,str_tim1r
     &                     ,radfact)
*--------              steered MD implemented (Jarzynsky on extra itors) 
                      call steer_dynamics(litor,nitors_added,strtors
     &                     ,potit(1,2),eqdied,eqdied1,tor_vel,rtime
     &                     ,str_time0,str_time1,str_tim0r,str_tim1r
     &                     ,radfact)
                      if(addstrcom) THEN 
                        if(rtime.gt.str_time0.and.rtime.le.str_time1)
     &                       THEN  
                          r0com_t=r0com+str_com*(rtime-str_time0)
                        end if
                      end if
                    else
c========           steering on an input define path in the space of
c========           bonds bends and tors (see read_potential) 
                      call steer_along_path(lstretch,nbonds_added,lbend
     &                     ,nbends_added,litor,nitors_added,strbonds
     &                     ,strbends,strtors,potbo(1,2),potbe(1,2)
     &                     ,potit(1,2),potbo(1,1),potbe(1,1)
     &                     ,potit(1,1),radfact,timerec,ntimes,pathbo
     &                     ,pathbe,pathto,activebo,activebe,activeto
     &                     ,rtime,eqdist,eqang,eqdied,str_vel,ben_vel
     &                     ,tor_vel,indexslice)
                    end if    
*-------            computes current values of steered coordinate 
                    call steer_variables(nbonds_added,nbends_added
     &                   ,nitors_added,atom_b1,atom_b2,atom_be1
     &                   ,atom_be2,atom_be3,atom_it1,atom_it2,atom_it3
     &                   ,atom_it4,xp0,yp0,zp0,bond_added,bend_added
     &                   ,tors_added,potit(1,2),litor,nitors_added)
                    if(addstrcom) THEN ! compute actual lig-tar COM distance 
                      laux=.false.
                      call com_restraint (kcom_md,r0com_t,xp0,yp0,zp0
     &                     ,mlig,mtar,mass,rcom,fpx_p,fpy_p,fpz_p,ucom
     &                     ,ilig1,ilig2,itar1,itar2,laux)
                    END IF
                  END IF
215               iaux=-100
#ifdef _MPI_ 
                  iaux=para_index
#endif
                  if(kplot_stn0.gt.0) THEN 
                    if(.not.path_steer.or.alchemy) THEN
                      call comp_work(nbonds_added,nbends_added
     &                     ,nitors_added,strbonds,strbends,strtors
     &                     ,bond_added,bend_added,tors_added,eqdist
     &                     ,eqdist1,eqang,eqang1,eqdied,eqdied1
     &                     ,force_cost,force_ang,force_died,str_vel
     &                     ,ben_vel,tor_vel,str_time0,str_time1
     &                     ,str_tim0r,str_tim1r,radfact,tn0,nplotstn0
     &                     ,ninn0,kplot_stn0,workbo,workbe,workto,rtime
     &                     ,iaux)
                      if(addstrcom) THEN 
                        call comp_work_com(rcom,r0com_t,kcom,str_com
     &                       ,rtime,str_time0,str_time1,tn0,workcom
     &                       ,diffcom,nplotstn0,ninn0,kplot_stn0,iaux) 
                      END IF  
                    else
                      if(indexslice(1).gt.0) THEN 
                        str_time0=timerec(indexslice(1))
                        str_time1=timerec(indexslice(2))
*-----                   no reverse path in steer_path
                        str_tim0r=0.0
                        str_tim1r=0.0
                        call comp_work(nbonds_added,nbends_added
     &                       ,nitors_added,strbonds,strbends,strtors
     &                       ,bond_added,bend_added,tors_added,eqdist
     &                       ,eqdist1,eqang,eqang1,eqdied,eqdied1
     &                       ,force_cost,force_ang,force_died,str_vel
     &                       ,ben_vel,tor_vel,str_time0,str_time1
     &                       ,str_tim0r,str_tim1r,radfact,tn0,nplotstn0
     &                       ,ninn0,kplot_stn0,workbo,workbe,workto
     &                       ,rtime,iaux)
                      end if
                    end if
                  end if
                  CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n0
     &                 ,fpy_n0,fpz_n0,ubond,uslvbon,ubend,uslvben,uitors
     &                 ,uslvitor,ndim)
                  ! compute force due to lig-tar COM restraint
                  if(addstrcom) THEN 
                    laux=.true.  
                    call com_restraint (kcom_md,r0com_t,xp0,yp0,zp0,mlig
     &                   ,mtar,mass,rcom,fpx_n0,fpy_n0,fpz_n0,ucom,ilig1
     &                   ,ilig2,itar1,itar2,laux)
                  end if
#ifdef _MPI_
                  if( rem_run .or. (sim_tempering .and.
     &                 .not.sim_tempering_sge) ) then
#else
                  if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
                    if(.not.rem_segment) THEN
                      call dscal(ntap,rem_factor(1),fpx_n0,1)
                      call dscal(ntap,rem_factor(1),fpy_n0,1)
                      call dscal(ntap,rem_factor(1),fpz_n0,1)
                    endif
                  endif

                  if(nadd_tpg.gt.0.or.addstrcom) THEN 
#ifdef  _OMP_
                    CALL comp_fcm_omp(nprot,ntap,iprot,fpx_n0,fpy_n0
     &                   ,fpz_n0,fcax_n0,fcay_n0,fcaz_n0,mass,tmass1,oc)
#else
                    CALL comp_fcm(nprot,protl,fpx_n0,fpy_n0,fpz_n0
     &                   ,fcax_n0,fcay_n0,fcaz_n0,mass,tmass,oc)
#endif                    
                  END IF

                  CALL correc(vpx,vpy,vpz,fpx_n0,fpy_n0,fpz_n0,mass,ntap
     &                 ,tn0)

                  IF(cnstpp .NE. 0) THEN
                    CALL rattle_correc(tn0,lx,ly,lz,vpx,vpy,vpz,ntap
     &                   ,cnstp,dssp,coeffp,cnstpp,mass,dnit,cnst_protp
     &                    ,cnst_protl,mim_lim,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
                  
*----------  COM correctionn due to addstrcom option   
                  
                  if(nadd_tpg.GT.0.or.addstrcom) THEN 
                    CALL correc(vcax,vcay,vcaz,fcax_n0,fcay_n0,fcaz_n0
     &                   ,tmass,nprot,tn0)
                    if(cpress) call comp_stress_add_tpg(nprot,stresst0
     &                   ,co,oc,xpcm,ypcm,zpcm,fcax_n0,fcay_n0,fcaz_n0
     &                   ,vco,masspp,tn0)
                  end if
                END DO
*---------   End of n0 Respa loop

                CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n1
     &              ,fpy_n1,fpz_n1,fudge,lj_fudge,abmd_dir,puhyd
     &              ,conf_bnd_slt_n1,conf_bnd_slv_n1,coul_bnd_slt_n1
     &              ,coul_bnd_slv_n1,unb14,cnb14,ungrp,cngrp,uptors
     &              ,uslvtor,mapdn,nmapdn,uumb,gr)
                dwrk14=dwrk
                if(alchemy.and.clewld) dwrk14t=dwrk14t+dwrk14
     &               /dfloat(n1respa)
#ifdef _MPI_
               if(rem_run .or. (sim_tempering .and. 
     &              .not.sim_tempering_sge)) then 
#else
               if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
                 if(.not.rem_segment) then  !if rem_segment is .TRUE. the n1 forces are scaled inside fnb14_rem and fptor_rem
                   call dscal(ntap,rem_factor(2),fpx_n1,1)
                   call dscal(ntap,rem_factor(2),fpy_n1,1)
                   call dscal(ntap,rem_factor(2),fpz_n1,1)
                 endif
               endif
               
*=======================================================================
*--- Call force routine for abmd on a torsion angle          -----------
*=======================================================================
               
               CALL correc(vpx,vpy,vpz,fpx_n1,fpy_n1,fpz_n1,mass,ntap
     &              ,tn1)

            END DO
            
*=======================================================================
*          |-->  END of n1-loop                                       ---
*=======================================================================
            
            IF(thermos) CALL verlet_free_eta(neta,eta,etap,tm)

            if(thermos.AND.steer_temperature)  call jartemp(t,temp,tiniz
     &           ,tfina,velt,rtime,timetiniz,timetfina,ucns_h,ucos_h
     &           ,ucns_l,ucos_l,ucns_m ,ucos_m,urcs_h ,urcs_l,urcs_m
     &           ,uslvbon,uslvben,uslvtor ,uslvitor ,conf_bnd_slv_n1
     &           ,coul_bnd_slv,coul_bnd_slv_n1 ,self_slv,efact,urcp_h
     &           ,urcp_l,urcp_m,ucop_h,ucop_l ,ucop_m,ucnp_h ,ucnp_l
     &           ,ucnp_m,ubond,ubend,uptors,uitors ,conf_bnd_slt_n1
     &           ,coul_bnd_slt,coul_bnd_slt_n1,self_slt ,ucek,pucek
     &           ,ucepr,uceh,cpress,urcsp_h,urcsp_l,urcsp_m ,ucosp_h
     &           ,ucosp_l,ucosp_m,ucnsp_h,ucnsp_l,ucnsp_m ,eer_m,eer_l
     &           ,eer_h,volume,pext,gascon,ninn0,nplottempstn0,tm
     &           ,kplot_steer_temp,nsteer_temp)


*=======================================================================
*---------- recomputes group position and scaled groups ---------------
*=======================================================================
                  
            CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
            
            CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga
     &           ,zpga)
            
*=======================================================================
*---        Computes M-forces at new coordinates                     ---
*=======================================================================
            
!           change the potential if alchemical transformations occur:
!           The lambda (added and removed) are changed according to the
!           protocol specified in the alchemical file (see ALCHEM dir) 
            if(time_steer.and.alchemy) THEN
!             find currect values of lambdas
              ladd0=ladd
              laddq0=laddq
              lrm0=lrm
              lrmq0=lrmq
              call alchemic(timerec,ntimes,ladded,laddedq,ladd,laddq
     &                ,lremoved,lremovedq,lrm,lrmq,rtime,indexslice)
!             set current values of atomic lambdas saving old values
!             to compute the work directly in mts_alchemy_forpp
              if(ladd.NE.1.D0.OR.laddq.NE.1.d0) THEN 
                do i=1,nat_added
                  lambda0(atom_added(i))=lambda(atom_added(i)) 
                  lambdaq0(atom_added(i))=lambdaq(atom_added(i)) 
                  lambda(atom_added(i))=ladd
                  lambdaq(atom_added(i))=laddq
                end do
                if(nat_added.ge.1) THEN 
                  dlam= lambda0(1)-lambda(1)
                  dlamq=lambdaq0(1)-lambdaq(1)
                END IF
              ENDIF
              if(lrm.NE.0.D0.OR.lrmq.NE.0.d0) THEN 
                do i=1,nat_removed
                  lambda0(atom_removed(i))=lambda(atom_removed(i))
                  lambda(atom_removed(i))=-lrm
                  lambdaq0(atom_removed(i))=lambdaq(atom_removed(i))
                  lambdaq(atom_removed(i))=-lrmq
                end do
                if(nat_removed.ge.1) THEN 
                  dlam=lambda(1)-lambda0(1)
                  dlamq=lambdaq(1)-lambdaq0(1)
                END IF
              ENDIF
              if(mod(ninner,nprint).eq.0) THEN 
                if(ladd.lt.1.d0.and.(.not.ladd_end).and.nat_added.gt.0)
     &               THEN
                  write(kprint,10987) (1.d0-ladd)*100.d0
10987             FORMAT(" LJ-Alchemical transformation (add) in", 
     &                 " progress:"
     &                 ,3x,f8.2,"% done..."  )
                ENDIF
                if(lrm.gt.0.d0.and.(.not.lrm_end).and.nat_removed.gt.0)
     &               THEN
                  write(kprint,11988) lrm*100.d0
11988             FORMAT(" LJ-Alchemical transformation (remove) in",
     &                 " progress:",f8.2,"% done..." ) 
                END IF
                
                if(laddq.lt.1.d0.and.(.not.laddq_end).and.nat_added.
     &               gt.0)THEN
                  write(kprint,10989) (1.d0-laddq)
     &                 *100.d0
10989             FORMAT(" QQ-Alchemical transformation (add) in", 
     &                 " progress:"
     &               ,3x,f8.2,"% done..."  )
                ENDIF
                if(lrmq.gt.0.d0.and.(.not.lrmq_end).and.nat_removed.
     &                 gt.0)THEN
                  write(kprint,11990) lrmq*100.d0
11990             FORMAT(" QQ-Alchemical transformation (remove) in",
     &               " progress:",f8.2,"% done..." ) 
                END IF
              END IF
            END IF
            
            CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga
     &           ,xpcma,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_m,ucos_m
     &           ,virs_m,virsp_m,ucnp_m,ucop_m,ucnsp_m,ucosp_m,fpx_m
     &           ,fpy_m,fpz_m,stressd_m,nnlpp1,nnlpp2)
            dwrkm=dwrk
            if(clewld.and.alchemy) THEN
              CALL ferrf_alchemy_add(ss_index
     &             ,alphal,chrge,lambdaq,lambdaq0,xp0,yp0,zp0,listqq
     &             ,nlistqq,falch_slt,aux,fpx_m,fpy_m,fpz_m,erf_corr
     &             ,erf_arr_corr,delew,rlew,dwrk_al,ntap,intramol)
              coul_bnd_slt=coul_bnd_slt_rl+falch_slt
            END IF
            if(time_steer.and.alchemy) wrk_alchemy=wrk_alchemy+dwrkm
     &           +dwrkl+dwrkh + dwrkrc_m+dwrkrc_l+dwrkrc_h + dwrk_al
     &           +dwrk14
!#ifdef DEBUGALCHEMY
c            write(106,10071) rtime,wrk_alchemy*1.d4,dwrkm*1.d4,dwrkl
c     &           *1.d4, dwrkh*1.d4, dwrkrc_m*1.d4,dwrkrc_l*1.d4,dwrkrc_h
c     &           *1.d4,dwrk_al*1.d4,dwrk14*1.d4,rshell
c10071       format (f10.2, 9G15.5,3x,A1)
!#endif
            IF(abmd .AND. (associate .OR. dissociate)) CALL
     &           comp_abmd_fdiss(abmd_dir,rspset,diss_list,spring,xpa
     &           ,ypa,zpa,co,fpx_m,fpy_m,fpz_m,uumb,gr)
            
            IF(clewld) THEN
              if(alchemy) THEN 
                CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma
     &               ,ypcma,zpcma,urcsp_m,urcs_m,urcp_m,virsp_m,virs_m
     &               ,virp_m,fpx_m,fpy_m,fpz_m,fsin14,gsin14,fsbend
     &               ,gsbend,fsbond,gsbond,fscnstr_slt,fscnstr_slv
     &               ,coul_bnd_slt_rl,coul_bnd_slv,rshell,rshk,eer_m
     &               ,stressr_m,fudgec,self_slt)
                dwrkrc_m=dwrk
              ELSE
                CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &               ,zpcma,urcsp_m,urcs_m,urcp_m,virsp_m,virs_m,virp_m
     &               ,fpx_m,fpy_m,fpz_m,fsin14,gsin14,fsbend,gsbend
     &               ,fsbond,gsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt
     &               ,coul_bnd_slv,rshell,rshk,eer_m,stressr_m,fudgec
     &               )
              END IF
            END IF
            
#ifdef _MPI_
            if(rem_run .or. (sim_tempering .and. 
     &           .not.sim_tempering_sge) ) then 
#else
            if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
              if(.not.rem_segment) then 
!if rem_segment is .TRUE. the real part of the nb forces are scaled
!inside mts_forpp_rem, while the reciprocal contribution is not scaled
                call dscal(ntap,rem_factor(3),fpx_m,1)
                call dscal(ntap,rem_factor(3),fpy_m,1)
                call dscal(ntap,rem_factor(3),fpz_m,1)
                do i = 1,3
                  do j = 1,3
                    stressd_m(i,j) = stressd_m(i,j)*rem_factor(3)
                    stressr_m(i,j) = stressr_m(i,j)*rem_factor(3)
                  enddo
                enddo
              endif
            endif
            
*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   m-shell forces --------------------------
*=======================================================================
            
#ifdef   _OMP_
            CALL comp_fcm_omp(nprot,ntap,iprot,fpx_m,fpy_m
     &                 ,fpz_m,fcax_m,fcay_m,fcaz_m,mass,tmass1,oc)
#else
            CALL comp_fcm(nprot,protl,fpx_m,fpy_m,fpz_m,fcax_m,fcay_m
     &           ,fcaz_m,mass,tmass,oc)
#endif     
       
*=======================================================================
*---        Scale velocities according to the thermostat velocities ----
*=======================================================================
            
*=========  Propagate exp(i L_x tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(cpress,slt_exist,slv_exist
     &              ,ss_point(1,1),ss_point(1,2),etap,tm4,nprot,vcax
     &              ,vcay,vcaz,vpx,vpy,vpz,vco)
            END IF
            
*========= Propagate exp(i L_s + L_u tm/2) -==================================

            CALL correc(vpx,vpy,vpz,fpx_m,fpy_m,fpz_m,mass,ntap,tm)

            CALL correc_stress(cpress,m_skin,nprot,co,oc,vco,vcax,vcay
     &           ,vcaz,fcax_m,fcay_m,fcaz_m,stressd_m,stressr_m,volume
     &           ,press_m,press_kin,pext,tmass,masspp,tm,tm2)

*=========  Propagate exp(i L_x tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(cpress,slt_exist,slv_exist
     &              ,ss_point(1,1),ss_point(1,2),etap,tm4,nprot,vcax
     &              ,vcay,vcaz,vpx,vpy,vpz,vco)
            END IF
                  
*---------- Propagate exp(i L_y tm/4) ---------------------------------
            IF(thermos) THEN
               IF(isostress) THEN
                 CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &                ,nboxcnst,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap
     &                 ,cnstp,dssp,coeffp,cnstpp,mass,dnit,cnst_protp
     &                 ,cnst_protl,mim_lim,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               CALL comp_thermos_forces(cpress,ndf_thermos,ntap,ss_index
     &              ,nprot,co,vpx,vpy,vpz,vcax,vcay,vcaz,vco,mass,tmass
     &              ,masspr,t,fth)
               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF

*---------- Propagate exp(i L_z tm/2) =================================
            IF(cpress .AND. coupl_mol) THEN
               CALL correc_matr(tm2,co,oc,vco,gmgp,vcax,vcay,vcaz,nprot)
            END IF
            

*---------- Propagate exp(i L_y tm/4) ---------------------------------
            IF(thermos) THEN
               CALL comp_thermos_forces(cpress,ndf_thermos,ntap,ss_index
     &              ,nprot,co,vpx,vpy,vpz,vcax,vcay,vcaz,vco,mass,tmass
     &              ,masspr,t,fth)
               CALL correc_etap(neta,etap,fth,qmass,tm4)
               CALL comp_thermos_energy(neta,ndf_thermos,t,qmass,eta
     &              ,etap,uceh,hpot,temph)
            END IF

*=======================================================================
*----------  Compute averages and do some analysis at time step M ------
*=======================================================================

            INCLUDE 'mtsmd_avg_inc.CPP.f'
         END DO
         
*=======================================================================
*---     <--|  END of M-loop                                        ----
*=======================================================================
         
*=======================================================================
*---     Computes L-forces at new coordinates                       ----
*=======================================================================
         
         rshell='l'
         CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &        ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_l,ucos_l,virs_l
     &        ,virsp_l,ucnp_l,ucop_l,ucnsp_l,ucosp_l,fpx_l,fpy_l,fpz_l
     &        ,stressd_l,nnlpp,nnlpp1)
         if(time_steer.and.alchemy) dwrkl=dwrk

         IF(clewld) THEN 
           if(alchemy) THEN
             CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &            ,zpcma,urcsp_l,urcs_l,urcp_l,virsp_l,virs_l,virp_l
     &            ,fpx_l,fpy_l,fpz_l,fsin14,gsin14,fsbend,gsbend,fsbond
     &            ,gsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt_rl
     &            ,coul_bnd_slv,rshell,rshk,eer_l,stressr_l,fudgec
     &            ,self_slt)
                dwrkrc_l=dwrk
           ELSE
             CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma
     &            ,urcsp_l,urcs_l,urcp_l,virsp_l,virs_l,virp_l,fpx_l
     &           ,fpy_l,fpz_l,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &           ,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv
     &            ,rshell,rshk,eer_l,stressr_l,fudgec)
           END IF
         END IF
#ifdef _MPI_
         if(rem_run .or. (sim_tempering .and. .not.sim_tempering_sge) ) 
     &        then
#else
         if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
           if(.not.rem_segment) then 
!if rem_segment is .TRUE. the real part of the nb forces are scaled
!inside mts_forpp_rem, while the reciprocal contribution is not scaled
             call dscal(ntap,rem_factor(3),fpx_l,1)
             call dscal(ntap,rem_factor(3),fpy_l,1)
             call dscal(ntap,rem_factor(3),fpz_l,1)
             do i = 1,3
               do j = 1,3
                 stressd_l(i,j) = stressd_l(i,j)*rem_factor(3)
                 stressr_l(i,j) = stressr_l(i,j)*rem_factor(3)
               enddo
             enddo
           endif
         endif
         
*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   l-shell forces --------------------------
*=======================================================================
         
#ifdef _OMP_
         CALL comp_fcm_omp(nprot,ntap,iprot,fpx_l,fpy_l
     &        ,fpz_l,fcax_l,fcay_l,fcaz_l,mass,tmass1,oc)
#else
         CALL comp_fcm(nprot,protl,fpx_l,fpy_l,fpz_l,fcax_l,fcay_l
     &        ,fcaz_l,mass,tmass,oc)
#endif
         
*---     corrects velocities
         
         CALL correc(vpx,vpy,vpz,fpx_l,fpy_l,fpz_l,mass,ntap,tl)
         
         CALL correc_stress(cpress,l_skin,nprot,co,oc,vco,vcax,vcay
     &        ,vcaz,fcax_l,fcay_l,fcaz_l,stressd_l,stressr_l,volume
     &        ,press_l,press_kin,pext,tmass,masspp,tl,tl2)

      END DO
      
*=======================================================================
*---  <--| END of L-loop                                            ----
*=======================================================================

      
*=======================================================================
*---     Computes H-forces (reciprocal +direct lattice)             ----
*---               at new coordinates                               ----
*=======================================================================
      
      rshell='h'
      CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_h,ucos_h,virs_h,virsp_h
     &     ,ucnp_h,ucop_h,ucnsp_h,ucosp_h,fpx_h,fpy_h,fpz_h,stressd_h
     &     ,nnlpp0,nnlpp)
      if(time_steer.and.alchemy) dwrkh=dwrk

      IF(clewld) THEN 
        if(alchemy) THEN 
          CALL mts_furier_alchemy(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma
     &         ,zpcma,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h
     &         ,fpy_h,fpz_h,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &         ,fscnstr_slt,fscnstr_slv,coul_bnd_slt_rl,coul_bnd_slv
     &         ,rshell,rshk,eer_h,stressr_h,fudgec,self_slt)
                dwrkrc_h=dwrk
        ELSE
          CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma
     &         ,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h,fpy_h
     &         ,fpz_h,fsin14,gsin14,fsbend,gsbend,fsbond,gsbond
     &         ,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell
     &         ,rshk,eer_h,stressr_h,fudgec)
        END IF
      END IF

      
c--   update also the m neighbor lists with a PHONY call to l mts_forces
      rcutl_save=rcutl
      rcutl=rcutm
      CALL mts_forces(nstep,'l',xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &     ,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &     ,stressd_p,nnlpp,nnlpp1)
      rcutl=rcutl_save


#ifdef _MPI_
      if( rem_run .or. (sim_tempering .and. .not.sim_tempering_sge) ) 
     &     then
#else
      if( sim_tempering .and. .not.sim_tempering_sge ) then
#endif
        if(.not.rem_segment) then 
!if rem_segment is .TRUE. the real part of the nb forces are scaled
!inside mts_forpp_rem, while the reciprocal contribution is not scaled
          call dscal(ntap,rem_factor(3),fpx_h,1)
          call dscal(ntap,rem_factor(3),fpy_h,1)
          call dscal(ntap,rem_factor(3),fpz_h,1)
          do i = 1,3
            do j = 1,3
              stressd_h(i,j) = stressd_h(i,j)*rem_factor(3)
              stressr_h(i,j) = stressr_h(i,j)*rem_factor(3)
            enddo
          enddo
        endif
      endif
      
*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   h-shell forces --------------------------
*=======================================================================
      
#ifdef  _OMP_
      CALL comp_fcm_omp(nprot,ntap,iprot,fpx_h,fpy_h,fpz_h,fcax_h
     &     ,fcay_h,fcaz_h,mass,tmass1,oc)
#else
      CALL comp_fcm(nprot,protl,fpx_h,fpy_h,fpz_h,fcax_h,fcay_h,fcaz_h
     &     ,mass,tmass,oc)
#endif      
*=======================================================================
*---  Corrects velocities                                            ---
*=======================================================================
      
      CALL correc(vpx,vpy,vpz,fpx_h,fpy_h,fpz_h,mass,ntap,time)

      CALL correc_stress(cpress,h_skin,nprot,co,oc,vco,vcax,vcay
     &     ,vcaz,fcax_h,fcay_h,fcaz_h,stressd_h,stressr_h,volume
     &     ,press_h,press_kin,pext,tmass,masspp,time,time2)
      
*=======================================================================
*---  Interaction with a stochastic bath                            ---
*=======================================================================
      
      IF(landersen) THEN
         CALL comp_vel_labframe(vpx,vpy,vpz,vpx,vpy,vpz,co,nprot,protl
     &        ,mass,tmass,vcax,vcay,vcaz)
         CALL collision(ntap,vpx,vpy,vpz,mass,nutime,t,time)
         CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &        ,vcaz)
      END IF


*=======================================================================
*----------  Dump restart file and do tests at timeste Hp --------------
*=======================================================================

      IF(debug) THEN
        DO i=1,ntap
          fpx(i)=fpx_n0(i)+fpx_n1(i)+fpx_m(i)+fpx_l(i)+fpx_h(i)
          fpy(i)=fpy_n0(i)+fpy_n1(i)+fpy_m(i)+fpy_l(i)+fpy_h(i)
          fpz(i)=fpz_n0(i)+fpz_n1(i)+fpz_m(i)+fpz_l(i)+fpz_h(i)
        END DO
        fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
        write(6,764) fstep,fpx_n0(1),fpx_n1(1),fpx_m(1),fpx_l(1)
     &       ,fpx_h(1)
764     FORMAT ('Forces on atom 1 =',f12.2,5e12.5)
      END IF
 
           
C--- invert velocities if time is equal to str_tim1r and steer is T.
      if(time_steer.and.abs(rtime-str_tim1r).lt.1.d-10) THEN 
c       invertv=.false. 
c        write(6, *) "VELOCITIES INVERTED!!!!!!!!"
c        do i=1,ntap
c          vpx(i)=-vpx(i)
c          vpy(i)=-vpy(i)
c          vpz(i)=-vpz(i)
c        end do
      end if
        
      if(bussi) THEN 
        CALL bussi_thermostat(vpx,vpy,vpz,vcax,vcay,vcaz,t,efact
     &       ,gascon,mass,tmass,massinfty,ebussi,friction,ntap,nprot,co
     &       ,cnstpp) 
        IF(cnstpp .NE. 0) THEN
          CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap,cnstp
     &         ,dssp,coeffp,cnstpp,mass,dnit,cnst_protp,cnst_protl
     &         ,mim_lim,iret,errmsg)
          IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
        END IF
!        CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
!     &       ,vcaz)
      END IF
      INCLUDE 'mtsmd_dump_inc.CPP.f'

*=======================================================================
*--------         MD LOOP ENDS   !!                    -----------------
*=======================================================================
      
#ifdef _MPI_
! replica exchange in REM
! if there is only one replica, no exchanges are attempted
      if(rem_run.and.ntrajectories.gt.1) then 
        para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
        call print_remdata(time,ninner,mrespa,lrespa,para0) 
        call rem_exchange(time,nrject,nstep,iproc,ntrajectories,t)
      endif
#endif

! SGE envirnment
      if ( sim_tempering ) then
        if(always_accept.and.go_save) THEN 
          CALL dcopy(9,co,1,co_s,1)
          CALL dcopy(9,oc,1,oc_s,1)
          CALL dcopy(ntap,xp0,1,xp0_s,1)
          CALL dcopy(ntap,yp0,1,yp0_s,1)
          CALL dcopy(ntap,zp0,1,zp0_s,1)
          CALL dcopy(ntap,lx,1,lx_s,1)
          CALL dcopy(ntap,ly,1,ly_s,1)
          CALL dcopy(ntap,lz,1,lz_s,1)
          CALL dcopy(ngrp,xpg,1,xpg_s,1)
          CALL dcopy(ngrp,ypg,1,ypg_s,1)
          CALL dcopy(ngrp,zpg,1,zpg_s,1)
          CALL dcopy(nprot,xpcma,1,xpcma_s,1)
          CALL dcopy(nprot,ypcma,1,ypcma_s,1)
          CALL dcopy(nprot,zpcma,1,zpcma_s,1)
          CALL dcopy(ntap,fpx_m,1,fpx_m_s,1)
          CALL dcopy(ntap,fpy_m,1,fpy_m_s,1)
          CALL dcopy(ntap,fpz_m,1,fpz_m_s,1)
          CALL dcopy(ntap,fpx_l,1,fpx_l_s,1)
          CALL dcopy(ntap,fpy_l,1,fpy_l_s,1)
          CALL dcopy(ntap,fpz_l,1,fpz_l_s,1)
          CALL dcopy(ntap,fpx_h,1,fpx_h_s,1)
          CALL dcopy(ntap,fpy_h,1,fpy_h_s,1)
          CALL dcopy(ntap,fpz_h,1,fpz_h_s,1)
          CALL dcopy(ntap,fpx_n1,1,fpx_n1_s,1)
          CALL dcopy(ntap,fpy_n1,1,fpy_n1_s,1)
          CALL dcopy(ntap,fpz_n1,1,fpz_n1_s,1)
          CALL dcopy(ntap,fpx_n0,1,fpx_n0_s,1)
          CALL dcopy(ntap,fpy_n0,1,fpy_n0_s,1)
          CALL dcopy(ntap,fpz_n0,1,fpz_n0_s,1)
          CALL dcopy(9,stressd,1,stressd_s,1)
          CALL dcopy(9,stressr,1,stressr_s,1)
          CALL dcopy(9,gmgp,1,gmgp_s,1)
          CALL dcopy(3,eta,1,eta_s,1)
          CALL dcopy(3,fth,1,fth_s,1)
          CALL dcopy(1,temp,1,temp_s,1)
          go_save=.false.
        END IF
        call sge_drive ( bond_added, bend_added, tors_added, radfact,
     &       nstep, iproc, uslvbon, uslvben, uslvitor, ubond, ubend,
     &       uitors, efact, uslvtor, conf_bnd_slv_n1, coul_bnd_slv_n1,
     &       uptors, conf_bnd_slt_n1, coul_bnd_slt_n1, ucns_h, ucos_h,
     &       ucns_l, ucos_l, ucns_m, ucos_m, urcs_h, urcs_l, urcs_m,
     &       coul_bnd_slv, self_slv, urcp_h, urcp_l, urcp_m, ucop_h,
     &       ucop_l, ucop_m, ucnp_h, ucnp_l, ucnp_m, coul_bnd_slt,
     &       self_slt, urcsp_h, urcsp_l, urcsp_m, ucosp_h, ucosp_l,
     &       ucosp_m, ucnsp_h, ucnsp_l, ucnsp_m, eer_m, eer_l, eer_h,
     &       hpot, pext, pucek, uceh, ucek, ucepr, volume, cpress )
        if(always_accept.and.use_saved) THEN 
          CALL dcopy(9,co_s,1,co,1)
          CALL dcopy(9,oc_S,1,oc,1)
          CALL dcopy(ntap,xp0_s,1,xp0,1)
          CALL dcopy(ntap,yp0_s,1,yp0,1)
          CALL dcopy(ntap,zp0_s,1,zp0,1)
          CALL dcopy(ntap,lx_s,1,lx,1)
          CALL dcopy(ntap,ly_s,1,ly,1)
          CALL dcopy(ntap,lz_s,1,lz,1)
          CALL dcopy(ngrp,xpg_s,1,xpg,1)
          CALL dcopy(ngrp,ypg_s,1,ypg,1)
          CALL dcopy(ngrp,zpg_s,1,zpg,1)
          CALL dcopy(nprot,xpcma_s,1,xpcma,1)
          CALL dcopy(nprot,ypcma_s,1,ypcma,1)
          CALL dcopy(nprot,zpcma_s,1,zpcma,1)
          CALL dcopy(ntap,fpx_m_s,1,fpx_m,1)
          CALL dcopy(ntap,fpy_m_s,1,fpy_m,1)
          CALL dcopy(ntap,fpz_m_s,1,fpz_m,1)
          CALL dcopy(ntap,fpx_l_s,1,fpx_l,1)
          CALL dcopy(ntap,fpy_l_s,1,fpy_l,1)
          CALL dcopy(ntap,fpz_l_s,1,fpz_l,1)
          CALL dcopy(ntap,fpx_h_s,1,fpx_h,1)
          CALL dcopy(ntap,fpy_h_s,1,fpy_h,1)
          CALL dcopy(ntap,fpz_h_s,1,fpz_h,1)
          CALL dcopy(ntap,fpx_n1_s,1,fpx_n1,1)
          CALL dcopy(ntap,fpy_n1_s,1,fpy_n1,1)
          CALL dcopy(ntap,fpz_n1_s,1,fpz_n1,1)
          CALL dcopy(ntap,fpx_n0_s,1,fpx_n0,1)
          CALL dcopy(ntap,fpy_n0_s,1,fpy_n0,1)
          CALL dcopy(ntap,fpz_n0_s,1,fpz_n0,1)
          CALL dcopy(9,stressd_s,1,stressd,1)
          CALL dcopy(9,stressr_s,1,stressr,1)
          CALL dcopy(9,gmgp_s,1,gmgp,1)
          CALL dcopy(3,eta_s,1,eta,1)
          CALL dcopy(3,fth_s,1,fth,1)
          CALL dcopy(1,temp_s,1,temp,1)
          CALL comp_fcm(nprot,protl,fpx_m,fpy_m,fpz_m,fcax_m,fcay_m
     &         ,fcaz_m,mass,tmass,oc)

          CALL comp_fcm(nprot,protl,fpx_l,fpy_l,fpz_l,fcax_l,fcay_l
     &         ,fcaz_l,mass,tmass,oc)

          CALL comp_fcm(nprot,protl,fpx_h,fpy_h,fpz_h,fcax_h,fcay_h
     &         ,fcaz_h,mass,tmass,oc)

          IF(nadd_tpg.gt.0.or.addstrcom) THEN
            CALL comp_fcm(nprot,protl,fpx_n0,fpy_n0,fpz_n0,fcax_n0
     &           ,fcay_n0,fcaz_n0,mass,tmass,oc)
          END IF
104       format(A20,f10.1,2f8.3,I10)
          use_saved=.false.
          CALL zero0(etap,neta)
          CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,linit
     &         ,massinfty)
          IF(cnstpp .NE. 0) THEN
            CALL rattle_correc(time,xp0,yp0,zp0,vpx,vpy
     &           ,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &           ,dnit,cnst_protp,cnst_protl,mim_lim,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
          END IF
          CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &         ,vcaz)
        END IF
      endif
      
*=======================================================================
*-------Add hills and print out metadynamics data ---------------------
*=======================================================================
      
      if(meta_run) then 
        if(nstep.gt.nrject) then
          if(mod(nstep*time,meta_ts).lt.1.d-10) then 
            call add_hill(nstep,time)
            ! results must be printed here 
          endif
        endif
      endif
      
      
      IF(nstep .LT. nrject) GOTO 100
      IF(nstep .LT. maxstp) GOTO 100
      

      IF(annealing) THEN
         grad_max=-1.0D0
         DO i=1,ntap
            fpx(i)=fpx_n0(i)+fpx_n1(i)+fpx_m(i)+fpx_l(i)+fpx_h(i)
            fpy(i)=fpy_n0(i)+fpy_n1(i)+fpy_m(i)+fpy_l(i)+fpy_h(i)
            fpz(i)=fpz_n0(i)+fpz_n1(i)+fpz_m(i)+fpz_l(i)+fpz_h(i)
            IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
            IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
            IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
         END DO
         write(kprint,20000) DABS(grad_max)
         CALL prtfrc(kprint,ngrp,grppt,nres,M1,prsymb,beta,xp0,yp0,zp0
     &        ,fpx,fpy,fpz)
      END IF

*=======================================================================
*     Write timing
*=======================================================================


      if(omp_timing)  THEN 
        if(iproc.eq.0.and.ntrajectories.gt.1) write(6,7001) time_in
     &       ,tcpu_in,time_cn,tcpu_cn,time_n0,tcpu_n0,time_n1,tcpu_n1
     &       ,time_m,tcpu_m,time_l,tcpu_l,time_h,tcpu_h,time_u,tcpu_u
     &       ,time_f,tcpu_f
        write(kprint,7001) time_in,tcpu_in,time_cn,tcpu_cn,time_n0
     &       ,tcpu_n0,time_n1,tcpu_n1,time_m,tcpu_m,time_l,tcpu_l,time_h
     &       ,tcpu_h,time_u,tcpu_u,time_f,tcpu_f
7001    FORMAT( 15x,23x,4x,"Elaps.",9x," CPU",/
     &          15x,"      time(integr)  = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(cnst)    = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(intra NO)= ",G15.5,G15.5," sec.",/   
     &          15x,"      time(intra N1)= ",G15.5,G15.5," sec.",/   
     &          15x,"      time(inter M) = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(inter L) = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(inter H) = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(updte U) = ",G15.5,G15.5," sec.",/   
     &          15x,"      time(recpr F) = ",G15.5,G15.5," sec.")   
      END IF
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      elaps= -elaps + elapse

      if(iproc.eq.0.and.ntrajectories.gt.1) THEN 
        write(6,*)
        write(6,60030)
        write(6,17000) gcpu
        write(6,18000) elaps
        write(6,60200) gcpu/time_fs
        write(6,60300) elaps/time_fs
        write(6,60030)
      END IF

      write(kprint,*)
      write(kprint,60030)
      write(kprint,17000) gcpu
      write(kprint,18000) elaps
      write(kprint,60200) gcpu/time_fs
      write(kprint,60300) elaps/time_fs
      write(kprint,60030)

      IF(wrtgyr) THEN
         CLOSE(kgyr)
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(////20('*'),'  M o l e c u l a r   T o p o l o g y  ',
     &     20('*')/////)
1100  FORMAT('*',13(' '),
     &     ' I n p u t   O p e r a t i o n s   C o m p l e t e d ',
     &     12(' '),'*'/'*',78(' '),'*')
1110  FORMAT('*',13(' '),
     &     '   M . D.   I n t e r m e d i a t e   R e s u l t s  ',
     &     12(' '),'*'/'*',78(' '),'*')
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
2000  FORMAT(/'<------ Linked cell list indexing updated -------->'/)
13000 FORMAT(/22x,' Temperature has been rescaled ',i5,' times '/)
14000 FORMAT(/22x,'Adjusting bond length to Force Field.'/
     &        22x,'     This will take a while...'//) 
20000 FORMAT(
     &     21x,'***************************************'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************')
70200 FORMAT(5x,'<------ Dumping Restart File ------->'/)
70300 FORMAT(//5x,'<------ Reading Restart File ------->',A)
70400 FORMAT(5x,'<------ Restart File Read in ------->'//)
70100 FORMAT('Velocities have been rescaled ---->'/)

17000 FORMAT(  15x,' Total cpu time for the run       = ',g15.5)
18000 FORMAT(  15x,' Total elapsed time for the run   = ',g15.5)
60200 FORMAT(  15x,' Averaged time per femtosecond    = ',3x,g15.5)
60300 FORMAT(  15x,' Averaged elapsed per femtosecond = ',3x,g15.5//)
15000 FORMAT(/// 10x,'* * * * r - R E S P A  i s  O N  * * * *'  //)
15500 FORMAT(/   10x,'======= timing ========================='  /)
16011 FORMAT(/5x,'CPU and elapse time for update          =',2G15.5)
15011 FORMAT(/5x,'CPU and elapse time for linked cell indexing ='
     &     ,2G15.5)
16001 FORMAT(/5x,'CPU  time for h-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5,
     &       /5x,'elps time for h-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5)
16002 FORMAT(/5x,'CPU  time for l-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5,
     &       /5x,'elps time for l-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5)
16003 FORMAT(/5x,'CPU  time for m-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5,
     &       /5x,'elps time for m-contribution: RECP =',G15.5,
     &     ' DIR =',G15.5,' TOT =',G15.5)
14004 FORMAT(/5x,'CPU  time for n1-contribution  =',G15.5, 
     &       /5x,'elps time for n1-contribution  =',G15.5) 
16004 FORMAT(/5x,'CPU  time for n0-contribution  =',G15.5,
     &       /5x,'elps time for n0-contribution  =',G15.5) 
10067 FORMAT(/5x,'THEORIC SPEED UP FOR NON BONDED PART =',f10.3)
10068 FORMAT(/5x,'OVERALL THEORIC SPEED UP =',f10.3/)
80033 format(/5x,'Expected ELAPSED time for the RUN:',I4,
     &     ' hours and ',I2, ' min',/,  
     &       /5x,' Expected average time per M step:',g15.5,' sec.'/  
     &       /5x,' Expected average time per femto :',g15.5,' sec.'/)  
60030 FORMAT(/10x,'==========================================='/)
78410 FORMAT
     &(//' *******ERROR: MAXT for PME too small. INCREASE MAXT.'//)
70700 FORMAT(//15x,' Program Stops Smoothly. Restart Dumped'//)
70120 FORMAT('Velocities of the barostat have been rescaled   ---->'/)
70130 FORMAT('Velocities of the thermostats have been rescaled   ---->'/
     &     )
20900 FORMAT(//' ****** WARNING ! WARNING ! WARNING ***************' /
     &        /' ******   drift_remove ON           ***************' /
     &        /' ****** WARNING ! WARNING ! WARNING ***************' //)
10072 FORMAT(//' --- Total energy remove per particle = ', f12.4,
     &         ' (KJ/mole)'/
     &         ' --- Number of dirty scaling          = ', I10/
     &         ' --- Frequency of scaling             = ',f10.3,
     &         ' 1/ps '/)
10977 FORMAT(//'*******WARNING: NO COFACTOR ATOMS SELECTED '/ 
     &     ' NCOFACTOR IS SET TO ZERO AND NO FIELD I COMPUTED'//)
80000 FORMAT('REMARK   Rigid body fit on CA atoms')
80100 FORMAT('REMARK   Rigid body fit on heavy atoms')
89000 FORMAT(//
     &     21x,'***************************************',/
     &     21x,'*                                     *',/
     &     21x,'*     Force Field Parameters for      *',/
     &     21x,'*       the current molecules         *',/
     &     21x,'*                                     *',/
     &     21x,'*                                     *',/
     &     21x,'***************************************'//)
89500 FORMAT(//21x,'***************************************'//)
90000 FORMAT('BOND')
91000 FORMAT('BENDINGS')
92000 FORMAT('TORSION PROPER')
93000 FORMAT('TORSION IMPROPER')
94000 FORMAT('END')
95000 FORMAT('#   Constrained stretchings follow')
96000 FORMAT('NONBONDED MIXRULE')
97000 FORMAT('#   Warning! Mixrules always assumed')
98000 FORMAT(/38('>'),38('<')/
     &     ' Warning number of groups in the restart file is ',
     &     'different from current'/
     &     5x,' Restart group number is ',i7,'  current is  ',i7/38('>')
     &     ,38('<')/)
99000 FORMAT(/38('>'),38('<')/
     &     ' Warning number of molecules in the restart file is ',
     &     'different from current.'/
     &     5x,' Restart molecule number is ',i7,'  current is
     &     ',i7/38('>'),38('<')/)
1011  FORMAT( //"# **** DATA FROM STEERED DYNAMICS *****", //
     &         '# The data in this file must be post-processed; '
     &        , 'Lines can be of the following type:  '/
     &        , ' '/
     &        , '#  a) bond n. I  TIME  R0(T) R  WORK(I) '/
     &        , '# where: I is the i-the steered bond; TIME is the '
     &        , ' actual time of the simulation; R0(T) is the driven'/
     &        , '# equilibrium coordinate; R is the actual bond coord' 
     &        , 'inate and WORK(I) is the work done on this coordinate'/ 
     &        , ' '/
     &        , '#  b) bend n. I  TIME  ALPHA0(T) ALPHA  WORK(I) '/
     &        , '# where: I is the i-the steered bend; TIME is the '
     &        ,'actual time of the simulation; ALPHA0(T) is the driven'/
     &        , '# equilibrium coordinate; ALPHA is the actual bend ' 
     &        , 'coordinate and WORK(I) is the work done on this '
     &        , 'coordinate'/ 
     &        , ' '/
     &        , '#  c) tors n. I  time  ALPHA0(T) ALPHA  WORK(I) '/
     &        , '# where: I is the i-the steered tors; TIME is the '
     &        ,'actual time of the simulation; ALPHA0(T) is the driven'/
     &        ,'# equilibrium coordinate; ALPHA is the actual bend cord' 
     &        ,'inate and work(i) is the work done on this coordinate'/ 
     &        ,' '/
     &        ,'#  d) WORKfunction T LAMBDA WBO WBE WTO TOTWORK '/
     &        ,' where T is the actual time of the simulation ; LAMBDA '
     &        ,' is the driven overall coordinate; WBO is the total' /
     &        ,'# work done on all the stretchings; WBE is the total '
     &        ,'work done on all the torsions; TOTWORK is the total '/
     &        ,'work done on the system during the nonequilibrium .' 
     &        ,'process. Units are: Kj/mol, angs and degree'//)

8000  FORMAT(//
     &'==============  R E P L I C A  E X C H A N G E  S I ',
     &'M U L A T I O N  ==========='// 
     &' ------ Current processor  is:    ',i10/
     &' ------ Number of batteries:      ',i10/
     &' ------ Number of processors:     ',i10/
     &' ------ Number of  replicas:      ',i10/
     &' ------ Lowest temperature:       ',f10.2/ 
     &' ------ Highest temperatures:     ',3f10.2/ 
     &' ------ Exchanges attempted each: ',f10.2,' fs',/ 
     &' ------ Print diagnostic each     ',f10.2,' fs',/
     &' ------ REM.set updating time     ',f10.2,' fs',/
     &' ------ Current replica is:       ',I10 / 
     &' ------ Actual  replica is:       ',I10 / 
     &' ------ Battery index is:         ',I10 / 
     &' ------ Current temperatures:     ',3f10.2 /
     &' ------ Skip force (no nb intramol):',L6 /
     &' ------ ', a40)

1003  FORMAT(//
     &'==============  FFTW LIBRARIES IS USED FOR PME',
     &'  ================================') 
1004  FORMAT(//
     &'==============  PUBLIC FFT LIBRARY IS USED  FOR PME',
     &'  ===========================') 

8001  FORMAT(/
     &' ------ Local scaling/solute tempering  is ON: definition',
     &     ' of  SOLUTE follows')

8002  FORMAT(
     &' ------ group n. ',i4, ' goes from atom ', i6, '(',A4,')',  
     &        ' to atom', i6, '(',A4,')')

8003  FORMAT(
     &' ------ Total number of atoms of the solute is ',i6)

8004  FORMAT(/' ------ Both solute-solvent and solute-solute are 
     & scaled',//
     x'====================================================',
     x'============================'//)

8005  FORMAT(/' ------ Only intrasolute interactions are scaled ',//
     x'====================================================',
     x'============================'//)

8006  FORMAT(/' ------ Only solute-solvent  interactions are scaled ',//
     x'====================================================',
     x'============================'//)

8007  FORMAT(/' =======   Moves are ALWAYS ACCEPTED:',/
     &       ,' -------   Target Replica is: ',i5, / 
     &       ,' -------   Cycle lenght Target Replica is: ',f10.1, //) 
1673  FORMAT(// ' ================= NON EQUILIBRIUM THERMAL CHANGES',
     &     ' ====================',/
     &     ' ======== Inizial temperature is: ' ,f10.3,/  
     &     ' ======== Final temperature is: ',f12.3,/
     &     ' ======== T-steer non equilibrium dynamics starts at:',f15.1
     &     ,' fs and ends at', f15.1,' fs') 
9000  FORMAT(//
     x'==============  M E T A D Y N A M I C S  S I ',
     x'M U L A T I O N  ==========='// 
     x' ------ ',A40/
     x' ------ Number of RCs            ',i10/ 
     x' ------ Hills deposition interval',1f10.4, ' fs'/ 
     x' ------ Hills height and width(s)',f10.4,5f10.4/) 

3434      FORMAT(/
     &         /" ====== Initial alchemic work due PBC reference state
     &:",f10.1,f12.4//)
3436      FORMAT(/
     &         /" ====== Final alchemic work due PBC reference state
     &:",f12.4//)
      RETURN
      END

!-----------------------------------------------------------------------------------
      subroutine check_coordinates(ntap,co,xpa,ypa,zpa,imin,jmin,laux
     &     ,drmin)
!-----------------------------------------------------------------------------------
      implicit none
      INCLUDE 'pbc.h'
! arguments
      integer kprint,ntap,imin,jmin
      real*8  xpa(*),ypa(*),zpa(*),drmin,co(3,3)
      logical laux   
! local 
      integer i,j
      real*8  dx,dy,dz,dr,xgg,ygg,zgg

!     executable statemnts
      drmin=1.d30
      imin=0
      jmin=0
      laux=.false. 
      do i=1,ntap-1   
        do j=i+1,ntap 
          xgg=(xpa(i)-xpa(j)) - 2.d0*PBC(xpa(i)-xpa(j))
          ygg=(ypa(i)-ypa(j)) - 2.d0*PBC(ypa(i)-ypa(j)) 
          zgg=(zpa(i)-zpa(j)) - 2.d0*PBC(zpa(i)-zpa(j)) 
          dx=co(1,1)*xgg+co(1,2)*ygg+co(1,3)*zgg
          dy=            co(2,2)*ygg+co(2,3)*zgg
          dz=                        co(3,3)*zgg
          dr = dx*dx+dy*dy+dz*dz
          if(dr.lt.drmin) THEN 
            drmin=dr 
            imin=i
            jmin=j
          END IF
          if (dr.LT.0.8d0) THEN
            laux=.true.
          END IF
        end do
      end do
      return
      end
