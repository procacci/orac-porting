module cpropar
  use parst
  implicit none
  INTEGER maxstp,nsave,nrject,nconf,nprint,nplot,nprop,nform,icl,icm &
       ,icn,nmol,nhskip,nupdte,nupdte_index,nunit,llacc,lldon,nbone &
       ,ntap,lbond,lbend,ltors,litor,int14p,nspw,maxrun,nflag(10),scale &
       ,kgyr,ngyr,int13p,dunitgi,dunfrci,dunitgr,dunfrcr,duitgi(10) &
       ,dufrci(10),ngrdon,ntop,lstretch,lstrtch(2,m9),lconstr,lcnstr(2,m9),nxrms,nvi &
       ,ndipole,nrespa,n1respa,n0respa,mrespa,lrespa,ntmtss(3),ntmtsp(5) &
       ,ncofactor,nsegs(2,m1),nfragm,fragm(2,m1),gofr_nprint &
       ,gofr_navg,gofr_ncomp,nmol_create,mim_lim,ndf_thermos(3),neta &
       ,navg_str,nrms,navg_str_xrms,atoms_diss(2,2),mol_diss(2) &
       ,diss_list(2,m1),start_anl,stop_anl,update_anl &
       ,top_bonds(ntopol+1),top_bendings(ntopol+1),top_ptors(ntopol &
       +1),top_itors(ntopol+1),ntop_print,nfragm_dist,nhbonds,nabmd &
       ,corr_atoms(1+m1),nbinew,nstep_freq,nbonds_added,atom_b1(500),atom_b2(500), &
       nbends_added,atom_be1(500),atom_be2(500),atom_be3(500),nitors_added,atom_it1(500) &
       ,atom_it2(500),atom_it3(500),atom_it4(500),nscalemass,ntimes,nbonds_steer_added & 
       ,nbends_steer_added,nitors_steer_added,nat_added,nat_removed,atom_added(1000) &
       ,atom_removed(1000),nlistqq,listqq(3,1000*1999),atomtoscale(m1)
  
  INTEGER protl(m1*2),nprot,protlb(m1*2),nprotb,nprot_charges,prot_charges(npm) &
       ,clsthl(m1*2),nclsth
  
  CHARACTER*1 shell_pme
  
  INTEGER lacc(2,m1),ldon(2,m1),nbtype(m1),nres(m1,2),nbun,mend(nores) &
       ,mback(m1),lbnd(2,m9),concta(m1,m10),concth(m1,m10) &
       ,lbndg(3,m2),ltor(4,m3),int14(2,m3),type14(m3),litr(4,m4),nhtype(m6,m6),lpnbd &
       ,lphyd,pea,peb,pec,nec,nascii,atres(2,nores),ngrp,grppt(2,m11),ingrp(2,m12),ingrpp &
       ,mres(2,nores),resg(m11),atomg(m1),int13(2,m2),bsitp(m13,m14),asitp(m13,m14),nbsitp(m14),nasitp(m14),nrigg,nrigat &
       ,wres,totvaf,prtvaf,novaf,divide_spline,iseed,protg(2*m11),groupp(m11) &
       ,atomp(m1),iz,nv1,nv2,nv3,nv4,nospt,ns1,ns2,ns3,ns4,niop,nios &
       ,ngprint,itor_ptype,nplot_fragm,nplot_steer,nplottempstn0 &
       ,nplotstn0,nplotalch0,nplot_center,ncx,ncy,ncz,type_table(m6,m6),abmd_tors_num &
       ,divide_records,atom_record
  INTEGER nplmap,plmap(40)
  INTEGER rgmap(m1),rgmapc(m1),nsct,rgnato,nfield,ffield(f1,f2),fexcl(f3,f2),nfatm(f2),nexatm(f2), &
       mfield,nesad(2),neead(2),maxsim,rggrup,rggmap(m11),menerg, &
       kenerg,nhoov,anntom,ancorr,annslv,annpro,anpoint(2,ndan1), &
       nvol,nfft1,nfft2,nfft3,pme_order,itprint,ifprint
  
  REAL*8  time,t,dtemp,radius,aaxis,baxis,caxis,volume, alphal,rkcut &
       ,rcut_corr,hrson,hrsoff,hanon,hanoff,rspon,mhoov(hoov),whoov &
       ,rspoff,rspcut,hrcut,hacut,wmtp,alf,bet,gam,phyd2(m5),phyd1(m5) &
       ,dtemph,delrg,gofr_cut,cut_field,alphaf,qmass(3),dtemppr,rcut_hb &
       ,acut_hb,a2cut_hb,annealing_fact,force_cost(500),eqdist(500),eqdist1(500) &
       ,force_ang(500),eqang(500),eqang1(500),force_died(500),eqdied(500),eqdied1(500),str_time0 &
       ,str_time1,str_tim1r,velscale(m1)
  REAL*8  rmol(3,8),co(3,3),oc(3,3),chrge(m1),mass(m1),potbe(m2,7), & ! maybe (m2,3)
       potto(m3,2),ptorj(m3,4),potit(m4,3),pnbd1(m6),pnbd2(m6) &
       ,pnbd3(m6),pnbd4(m6),ecc6(m5),ecc12(m5),ecc146(m5),ecc1412(m5),c6pw(m7),c12pw(m7),deps(m7),potbo(m9,2) &
       ,potbo_cnst(m9,2),ett,emass,echrg,ecut,xecm,yecm,zecm,lmbdt,tlmbdb &
       ,csecp(4,m13,m14),stret,epcut,omega,rspset,spring,pmass(m1),uealfa,uemin &
       ,uemax,xet(2),yet(2),zet(2),sjorg(m6),ejorg(m6),c6jorg(m5),c12jorg(m5),fudge,lj_fudge &
       ,lj_fudgeb,antstar,anfact,radii(m5),eps_energy,xpt0(m1),ypt0(m1),zpt0(m1) &
       ,erfc_bin,erfc_arr(4,mspline),rpass,nutime,erf_arr_corr(4,mspline),delew,rlew,ruew &
       ,hstep_freq,tiniz,tfina,timetiniz,timetfina,dock_radius,vsite
  REAl*8  pathbo(100,500),pathbe(100,500),pathto(100,500),timerec(100) &
       ,ladded(100),lremoved(100),laddedq(100),lremovedq(100),kcom,r0com,r1com
  INTEGER  activebo(100,500),activebe(100,500),activeto(100,500),ss_alchemy(m1)
  REAL*8  densurfpt,probrad
  REAL*8  xsct(10),ysct(10),zsct(10),rdsct(10),hdmass,wcoef,duitgr(10),dufrcr(10),massinfty
  REAL*8  pext,masspr,wpr
  REAL*8  taut,taup,mu,compressibility
  REAL*8  rcuth,rcutl,rcutm,rtolh,rtoll,rtolm,rneil,rneim,rneih,kl &
       ,klt,km,kmt,distmax,abmd_tors_dir,abmd_cryst_dir &
       ,abmd_kvect,lrm,ladd,lrmq,laddq,lrm0,ladd0,lrmq0,laddq0,lambda(m1),lambda0(m1),lambdaq(m1),lambdaq0(m1) &
       ,rmin(m1),epsm(m1),dhlam,dhlamq,wrk_alchemy,dwrk,friction &
       ,qshift,ljshift,qsoft,ljsoft
  LOGICAL linser,protei,solven,hydbnd,cutoff,mixrww,mixrwp,clewld, &
       stinit,debug,debug_rs,debug_ct,debug_st,debug_bt,debug_it &
       ,debug_pt,tpgfil,prttpg,electr,elinit,sction,hmass,tpgwbn &
       ,grpcut,rescm,lenerg,uenerg,amphi,bending,slvpdb,stoprun,pme &
       ,erfc_spline,erfc_spline_corr,ltest_times,md_respa &
       ,linked_cell,avg_str,abmd,abmd_tors,abmd_cryst,abmd_unbias &
       ,dissociate,diss_mol,diss_atoms,associate,remove_momentum &
       ,avg_rms,diffusion,avg_ca,avg_he,write_ff_pars,annealing &
       ,change_cell,erf_corr,write_grad,lform_field,scalemass &
       ,path_steer,always_scale_nose,alchemy,linitq,lrmq_end,print_qwork,do_workq,read_remark,intramol,ljintra
  LOGICAL lfield,prtatl,prthbl,prtnbl,prtcnl,prtbal,prtptl,prtitl &
       ,prtbndl,prtseq,prtipw,prtn14,fold,scharge,nonbnd,noneq &
       ,wrtgyr,iobdmp,dmprnd,hoover,thermos,widom,duitgl,dufrcl,vacf &
       ,time_corr,not_time_corr,analys,anxrms,anxrms_cell,anrms &
       ,anslv,anxca,anxbc,anxhe,anxal,anprot,adihed,stretch,cpress &
       ,pressure,isostress,isostressxy,coupl_mol,coupl_grp,coupl_atm,lberendsen &
       ,protvold,protvolc,pvexclud,sgroup,recstrc,mdsim,minimize &
       ,frequencies,conj_grad,steepest,igr,prttopl,prttor,prtben &
       ,prtstr,stretch_heavy,lslvint,template,spectra,clean &
       ,start_conf,energy_then_die,lgr,vacf_mts,replicate,gofr_avg &
       ,gofr_neighbor,gofr_intra,gofr,restart_cont,read_co,landersen &
       ,adjust_cnstr,tpgprm_read,occupy_space,fragm_dist,hbonds_res,hbonds_tot &
       ,ascii_nocell,xyz_nocell,dynamic &
       ,print_unformatted &
       ,steer_temperature,steer_opt,save_all_restarts,agbnp,lsasa & 
       ,check_coord,l_mpi,indtar(m1),indlig(m1),addstrcom,steer_com,bussi,moveitors &
       ,docking,gen_itp,shifted
  INTEGER  strbonds(500),strbends(500),strtors(500),ilig1,ilig2,itar1,itar2
  LOGICAL remltor(m3),remitor(m4),remint14(m3),rematom(m1)
  CHARACTER*80 titld0,gerald,dmpfil,cgroup,restart_file,restart_out,pdbfile
  CHARACTER*7 beta(m1),betb(m1),alnbd(m6)
  CHARACTER*8 anatom(ndan1)
  CHARACTER*8 prsymb(n5) 
  CHARACTER*80 title0
  INTEGER nbun_slv,nato_slv,ngrp_slv,lbond_slv,lbend_slv,ltors_slv &
       ,litor_slv,nbone_slv,int14p_slv &
       ,int13p_slv,llacc_slv,lldon_slv,mback_slv(slvatm),nres_slv(slvatm,2),mend_slv(slvunit) &
       ,grppt_slv(2,slv_group),int14_slv(2,slv4),type14_slv(slv4),int13_slv(2,slv3) &
       ,lacc_slv(2,slvatm),ldon_slv(slvatm) &
       ,nhtype_slv(types,types),lbnd_slv(2,slv2),lbndg_slv(3,slv3),ltor_slv(4,slv4) &
       ,litr_slv(4,slv5),atres_slv(2,slvunit) &
       ,concta_slv(slvatm,slv6),nbtype_slv(slvatm),bsitp_slv(slv_cns1,slv_cns2) &
       ,asitp_slv(slv_cns1,slv_cns2),nbsitp_slv(slv_cns2) &
       ,nasitp_slv(slv_cns2),nrigg_slv,nrigat_slv,icl_slv,icm_slv,icn_slv &
       ,ss_point(m1+1,2),ss_index(m1),ifixed(m1),nfixed,ndocking,jligand
  REAL*8  chrge_slv(slvatm),mass_slv(slvatm),potbo_slv(slv2,2),potbe_slv(slv3,7), & ! maybe (slv3,3)
       potto_slv(slv4,4),potit_slv(slv5,4) 
  
  CHARACTER*7 beta_slv(slvatm),betb_slv(slvatm),alnbd_slv(slvatm)
  CHARACTER*8 slv_residue,slv_redef_type
  LOGICAL slv_exist,slt_exist,slv_generate,slv_create,slv_redef &
       ,slv_read,adihed_slv,pdb_read,slt_create,slv_randomize &
       ,slv_add,lskipm,lskipl,lskiph

end module cpropar
