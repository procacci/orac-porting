      SUBROUTINE run_minimize(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm
     &     ,ypcm,zpcm)

************************************************************************
*                                                                      *
*     MTSMD is the driver of the MD run when multiple time scales      *
*     methods are used. The integration algorithm is r-RESPA           *
*     of order o(dt^3) (Tuckermann et al. JCP 97 1990 (1992))          *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi                                        *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*   RUN_MINIMIZE externals:    	       	       	       	       	       *  
*       add_energies appbou change_frame chkewl comp_dynamic_mat       *
*       comp_molmass comp_vel_labframe copy_protl dumprs	       *
*       erfc_spline_init erf_corr_cutoff fft_pme_init fndgrp	       *
*       get_total_energy inicmp int_corr_erf_spline *
*       lc_index lc_list linmin_total matinv mts_forces		       *
*       plotc prtacc prtat prtba prtbnd prtcn			       *
*       prtfrc prtit prtite_min prtpt prtsq readrs		       *
*       timer tr_inbox write_pot_bond write_pot_nbond xerror	       *
*       zero zero3x3 zeroa					       *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar
      use spme
      use sasamod

#ifdef _OMP_
      use omp_integr
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
      INTEGER mapnl(m8t,nthr)
#else
      INTEGER mapnl(*)
#endif
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),eta(*),xpcm(*)
     &     ,ypcm(*),zpcm(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'lc_list.h'
      
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif

*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      REAL*8   ddot,dnrm2
      LOGICAL near0,ok
      EXTERNAL  near0,ddot,dnrm2,f1dim_der

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      character*1  rshk
      
      INTEGER    kp,kt,nato_slt,iter_avg,idock,atom_prot(numpr),jlig
     &     ,itry,nstart,nend,nn,nlig

      INTEGER i,j,nsstt,tstep,nstep,iret
      INTEGER flag
      LOGICAL lupdate,nosave,lsasarun
      REAL*8  ucns,ucos,ucnsp,ucosp,ucnp,ucop,fudgec,fstep
      REAl*8 urcs,urcp,urcsp,eer,stressr(3,3),stressd(3,3),stressd_p(3,3
     &     ),stress_conf(3,3),stress_tot(3,3),prt(3,3)
     &     ,press_conf,press_kin,errca,errhe,errbc
     &     ,erral,drpca,drpbc,drphe,drpal,massinfty0
      real*8  xpcmr(npm),ypcmr(npm),zpcmr(npm)

      REAL*8  temp,elapse,puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors
     &     ,uconf,ucoul,ureal,urecp,ucek,pucek,tempt,fsrtal,tempr,temppr
     &     ,tcm,rcm,elaps,rms,unb14,cnb14,fsbond,fsbend,fsin14,purecp
     &     ,vfcp,tfcp,uumb,upconf,upcoul,gcpu,temph,uceh,hpot
     &     ,ubond,aux,hstep,uslvbon,uslvben,ucepr,gcpu_u,enisolv,vol_lig
      REAL*8  gr,pueng
      INTEGER navg
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
     &     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor,ssm_pnbd
     &     ,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul,ssm_epconf
     &     ,ssm_presst,ssm_press,ssm_pressc,ssm_pressk,ssm_co(3,3)
     &     ,ssm_st(3,3),ssm_volume,ssm_pv,ssm_temppra,ssm_eslvint
     &     ,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors

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

      REAL*8  temppra,energy

*----------- ARRAYS USED BY COMP_DYNAMIC_MAT ---------------------------
      
      REAL*8 fpx2(m1,cheb_order),fpy2(m1,cheb_order),fpz2(m1,cheb_order)
     &     ,d_mat(n_mat*(n_mat+1)/2),wk(n_mat*(n_mat+1)/2+n_mat)
     &     ,eigvl(n_mat),eigvc(n_mat,n_mat)
      INTEGER mad,mbd
      
*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER ngrp_old,nprot_old,nind(2) 
      INTEGER indxi(2,indmax),indxj(2,indmax),indxk(2,indmax)

      INTEGER, allocatable :: nnlpp0(:,:),nnlpp(:,:),nnlpp1(:,:)

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
      REAL*8  gpx(mb),gpy(mb),gpz(mb),hpx(mb),hpy(mb),hpz(mb),fpx(mb)
     &     ,fpy(mb),fpz(mb),fpx1(mb),fpy1(mb),fpz1(mb),vpx(mb),vpy(mb)
     &     ,vpz(mb),vpx1(mb),vpy1(mb),vpz1(mb),vcax(mb),vcay(mb),vcaz(mb
     &     )
      REAL*8  etap(hoov),vh1(hoov),vco(3,3)
      REAL*8  zz1(3,3),lzz,lzz1
      REAL*8  xpo(mb),ypo(mb),zpo(mb),co2(3,3),oc2(3,3)
      REAL*8  xp1(mb),yp1(mb),zp1(mb),xpa(mb),ypa(mb),zpa(mb),xpga(m11)
     &     ,ypga(m11),zpga(m11),xpcma(npm),ypcma(npm),zpcma(npm)
     &     ,tmass(numpr),tmassb(numpr),epslj(mb),sigm(mb)
      REAL*8  xp0_s(mb),yp0_s(mb),zp0_s(mb),ranf,dummy
      external ranf
      INTEGER dthree
      INTEGER numatoms,igrn,krdf(maxint*g1)
      INTEGER mapdn(2,mf),nmapdn(mf),tipo(mb),ih(mb),idum,jdum

*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external


*     Phony forces and energies for neighbor list

      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p(1),fpy_p(1),fpz_p(1),conf_bnd_slt,coul_bnd_slt
     &     ,conf_bnd_slv,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor
     &     ,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1
      
      LOGICAL dpress,dhoover,grflag,laux
      REAL*8  grad_max
      REAL*8  work(mspline),utotal,gg,dgg,gamma,eps2,eps,fret,utotal_old
     &     ,dutotal,utotalmin,useparated,ucomplex,ucomplext,enerdiss,tds
     &     ,volume_term

      INTEGER count

*==================== EXECUTABLE STATEMENTS ============================

      nosave=.false.
      lsasarun=.false.

*===  Check if the dimension of the work array are sufficient 


      IF(mb.LT.ntap) THEN
         errmsg=' While in MTSMD: PARAMETER MB dimensions the work'
     &        //' arrays is sufficient. Abort. '
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF

*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

*===  set few variable to zero


      flag=nflag(1)
      hstep=0.3D0
      urcs=0.0D0
      urcp=0.0D0
      urcsp=0.0D0
      eer=0.0D0
      rcutl=0.0D0
      rtoll=0.0D0
      ucek=0.0D0
      pucek=0.0D0
      uceh=0.0D0
      ucepr=0.0D0
      t=0.0D0
      pueng=0.d0
      puhyd=0.d0
      CALL zero3x3(stressr)
      CALL zero3x3(stressd)
      CALL zeroa(fpx,fpy,fpz,ntap,1)


      uumb=0.0D0
      conf_bnd_slv=0.0D0
      conf_bnd_slt=0.0D0
      tstep=0
      fudgec=1.d0-fudge
      nsstt=0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv
      iter_avg=0
      eps=eps_energy
      eps2=eps**2

#ifdef _OMP_
      mpp8=1+mpp/nthr
      allocate(nnlppf(mpp8,nthr),stat=iret)
      allocate(nnlpp0(mpp8,nthr),stat=iret)
      allocate(nnlpp(mpp8,nthr),stat=iret)
      allocate(nnlpp1(mpp8,nthr),stat=iret)
#else
      allocate(nnlpp0(mpp,1),stat=iret)
      allocate(nnlpp(mpp,1),stat=iret)
      allocate(nnlpp1(mpp,1),stat=iret)
#endif

      CALL comp_molmass(nprot,protl,mass,tmass)

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
         WRITE(kprint,70300)
         CALL readrs(kdump,nstep,temp,ntap,ngrp_old,nprot_old,xp0,yp0
     &        ,zp0,vpx,vpy,vpz,xpg,ypg,zpg,xpcma,ypcma,zpcma,vpx1,vpy1
     &        ,vpz1,vcax,vcay,vcaz,eta,etap,vh1,dthree,dhoover,sumarray
     &        ,ssmarray,navg,navg,co,oc,dpress,vco,zz1,lzz,lzz1,grflag
     &        ,krdf,igrn,maxint*g1,restart_cont)
         WRITE(kprint,70400)

         temp=0.0D0
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

      WRITE(kprint,'(//)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1100)
      WRITE(kprint,1100)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(//)')

      CALL zero(sumarray,ssmarray,navg,navg)
      nstep=0

*=======================================================================
*----- Save starting coordinates if docking=.true.
*=======================================================================

      if(docking) THEN 
        if(agbnp) THEN
          errmsg=
     &         " Docking only with SASA; change AGBNP to SASA in main in
     &put  "  
          CALL xerror(errmsg,80,1,2)
          stop
        endif
          
        utotalmin=1.D30
        ucomplex=1D30
        idock=0
        itry=0
        nprint=10000000
        call assign_sasa_types(ntap,nbtype,betb,kprint,sasa_nbtype
     &     ,sasa_betb,rvdw,psasa,ssasa,r2cut,rsolv)
        CALL dcopy(ntap,xp0,1,xp0_s,1)
        CALL dcopy(ntap,yp0,1,yp0_s,1)
        CALL dcopy(ntap,zp0,1,zp0_s,1)
        call find_ligand(protl,nprot,atom_prot,jlig)
        nn=0
!       reset ligand index to input value
        if(jligand.gt.0) THEN
          jlig=jligand 
        END IF
        do j=1,nprot
          if(j.eq.jlig) then 
            nstart=nn+1
            nend=nn+atom_prot(j)
          end if
          nn=nn+atom_prot(j)
        end do
        nlig=atom_prot(jlig)
        call volume_lig(nstart,nend,rvdw,sasa_nbtype,xp0_s,yp0_s,zp0_s
     &       ,vol_lig,nlig)
        write(kprint,1005) jlig,nlig,vol_lig,nprot,dock_radius
     &       ,ndocking
!         check coordinates before starting 
        if(idock.gt.0) write(kprint,10980) idock
10980   format(5x,/// " ====> Docking started round ---> ", i5)

10999   if(itry.gt.0.or.idock.eq.0) call rdocking(ss_index,xp0_s,yp0_s
     &       ,zp0_s,xpcm,ypcm,zpcm,mass,nprot,protl,nstart,nend,jlig
     &       ,nlig,ntap,dock_radius,xp0,yp0,zp0,idock) 


        CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
        call check_coordinates(ntap,co,xpa,ypa,zpa,idum,jdum
     &       ,laux,dummy)
        if(idock.gt.0) write(kprint,10789) dummy,idum,jdum
        if(laux.and.itry.lt.1000) THEN 
          itry=itry+1
          if(idock.gt.0) write(kprint,*)
     &         ' ====> Drug Coordinates are regenerated ',
     &         " -- tentative ", itry
          go to 10999 
        end if
        CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &           ,nprot,protl)
        CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
        CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &       ,ypcma,zpcma)
        CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
        CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo
     &       ,zpo)
        fstep=0.d0
        CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &       ,idum)
        call force_sasa(ntap,sasa_nbtype,sasa_betb,xpa,ypa,zpa,co,r2cut
     &       ,rsolv,rvdw,ssasa,psasa,sasaener,fpx,fpy,fpz,sasa) 
      end if

*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot,protl)

*=======================================================================
*---- Calculate group position  ----------------------------------------
*=======================================================================

      CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
      
*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================

      CALL fndgrp(nprot,protl,atomp)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

      IF(nflag(1) .NE. 0 .AND. ngrp_old .NE. ngrp ) THEN
         WRITE(kprint,98000) ngrp_old,ngrp
      END IF
      IF(nflag(1) .NE. 0 .AND. nprot_old .NE. nprot ) THEN
         WRITE(kprint,99000) nprot_old,nprot
      END IF

*=======================================================================
*----- Change the cell according to input commands ---------------------
*=======================================================================

      IF(change_cell) THEN
         CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &        ,zpcma)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)

         CALL dcopy(9,co2,1,co,1)
         CALL dcopy(9,oc2,1,oc,1)

         CALL change_origin(-1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_frame(co,oc,1,nprot,xpcma,ypcma,zpcma,xpcm,ypcm
     &        ,zpcm)
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
            CALL fft_pme_init(numatoms,nfft1,nfft2,nfft3,pme_order
     &           ,bsp_mod1,bsp_mod2,bsp_mod3,planf,planr)
#else
            CALL fft_pme_init(numatoms,nfft1,nfft2,nfft3,pme_order
     &           ,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
     &           ,bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork)
            if ( siz_Q .GT. MAXT ) THEN
               write(kprint,78410)
               stop
            END IF
#endif
            rshk=shell_pme
            IF(erf_corr) THEN
               CALL erf_corr_cutoff(oc,delew,rkcut,nfft1,nfft2,nfft3)
               CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &           ,erf_arr_corr,work)
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
     &           ,work,iret,errmsg,rkcut,erfc_spline_corr)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF
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
          nnlppf(1+j,i)=0
        end do
      end do
#else      
      nnlppf(1)=ngrp
      DO j=1,ngrp
         nnlppf(1+j)=j
      end do
#endif

*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================


      CALL prtite_min
      IF(prttpg) THEN
         WRITE(kprint,1000)
         WRITE(kprint,'(////)')
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
         WRITE(kprint,'(////77(''=''))')
      END IF

*======================================================================
*==== Write force field when required. Only the parameters used -------
*==== in the actual simulation will be printed ------------------------
*======================================================================

      IF(write_ff_pars) THEN
         WRITE(kprint,89000)
         write(kprint,*) "--> Begin PRM"
         WRITE(kprint,90000)
         CALL write_pot_bond(lstretch,2,lstrtch,m9,2,potbo,betb,1)
         WRITE(kprint,95000)
         CALL write_pot_bond(lconstr,2,lcnstr,m9,2,potbo_cnst,betb,-1)
         WRITE(kprint,94000)
         WRITE(kprint,91000)
         CALL write_pot_bond(lbend,3,lbndg,m2,4,potbe,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,92000)
         CALL write_pot_bond(ltors,4,ltor,m3,2,potto,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,93000)
         CALL write_pot_bond(litor,4,litr,m4,3,potit,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,96000)
         WRITE(kprint,97000)
         CALL write_pot_nbond(ecc6,ecc12,ecc146,ecc1412,nbtype,mass
     &        ,lj_fudgeb,betb,ntap)
         WRITE(kprint,94000)
         write(kprint,*) "--> End PRM"
         WRITE(kprint,89500)
      END IF


      if(agbnp) THEN
        call assign_agbnp_types(ntap,nbtype,epslj,sigm,alnbd
     &     ,pnbd1,pnbd2,tipo,kprint)
         do i=1,ntap
           ih(i) = 0
        end do
        call agbnpf_init(ntap,xp0,yp0,zp0,tipo,chrge,epslj,sigm,ih,iret)
      else if (lsasa.and.(.not.docking)) THEN 
        call assign_sasa_types(ntap,nbtype,betb,kprint,sasa_nbtype
     &       ,sasa_betb,rvdw,psasa,ssasa,r2cut,rsolv)
      end if
 
*======================================================================
*==== Write a Banner Before intermediate results ----------------------
*======================================================================

      WRITE(kprint,'(//)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1110)
      WRITE(kprint,1110)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(//)')

      if(gen_itp.and.nprot.eq.1.and.nmol.eq.0) THEN 
        call gen_itp_file
      end if
      IF(rneih.GT.0.d0) THEN
         lupdate=.true.
      else
         nnlpp0(1,1)=ngrp
         DO j=1,ngrp
            nnlpp0(1+j,1)=j
         end do
         lupdate=.false.
      END IF   

*======================================================================
*==== Add extra bonds (option ADD_BONDS in &POTENTIAL) 
*======================================================================

      if(nbonds_added.gt.0) THEN 
        WRITE(kprint,1210) nbonds_added
1210    FORMAT (/ ' --  Added ',i4, ' extra bonds (MTSMD LEVEL) '/
     &            13x,4x,'atom1',4x,'atom2',5x,'K',8x,'D0')
        do i=1,nbonds_added
          if(atom_b1(i).lt.nato_slt) THEN
            lstrtch(1,i+lstretch)=atom_b1(i)
          ELSE
            errmsg = 'ADD_BONDS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          END IF
          if(atom_b2(i).lt.nato_slt) THEN
            lstrtch(2,i+lstretch)=atom_b2(i)
          ELSE
            errmsg = 'ADD_BONDS atom out of solute range'
            CALL xerror(errmsg,80,1,1)
            STOP
          ENDIF
          potbo(i+lstretch,1) =  1000.0D0*force_cost(i)*4.184/(unite
     &         *avogad)
          potbo(i+lstretch,2) =  eqdist(i)
          WRITE(kprint,1201) lstrtch(1,i+lstretch),lstrtch(2,i+lstretch)
     &         ,force_cost(i),eqdist(i)
1201      FORMAT(10x,2i10,2f10.3)  
        end do
        lstretch=lstretch+nbonds_added
      end if

      if(nbends_added.gt.0) THEN 
        WRITE(kprint,1310) nbends_added
1310    FORMAT (/ ' --  Added ',i4, ' extra bendings (MTSMD LEVEL) '/
     &            13x,4x,'atom1',4x,'atom2',4x,'atom3',5x,'K',8x,'A0')
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
          pi = dacos(-1.d0) 
          potbe(i+lbend,1) = 1000.0D0*force_ang(i)*4.184/(unite
     &         *avogad)
          potbe(i+lbend,2) =  eqang(i)*(pi/180.d0)
          WRITE(kprint,1301) lbndg(1,i+lbend),lbndg(2,i+lbend)
     &         ,lbndg(3,i+lbend),force_ang(i),eqang(i)
1301      FORMAT(10x,3i10,2f10.3)
        end do
        lbend=lbend+nbends_added
      end if

c---  add extra itors given in READ_POTENTIAL using ADD_TORSIONS cmd
 
      if(nitors_added.gt.0) THEN 
        WRITE(kprint,1410) nitors_added
1410    FORMAT (/ ' --  Added ',i4, ' extra bendings (MTSMD LEVEL) '/
     &       13x,4x,'atom1',4x,'atom2',4x,'atom3',4x,'atom4',5x,'K',8x
     &       ,'A0')
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
          pi = dacos(-1.d0) 
          potit(i+litor,2)=eqdied(i)*(pi/180.0D0)
          potit(i+litor,1)=1000.0D0*force_died(i)*4.184/(unite*avogad)
          potit(i+litor,3) = 1.d0
          WRITE(kprint,1401) litr(1,i+litor),litr(2,i+litor) ,litr(3,i
     &         +litor),litr(4,i+litor),force_died(i),eqdied(i)
1401      FORMAT(10x,4i10,2f10.3)
        end do
        litor=litor+nitors_added
      end if

*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

1001  CONTINUE    ! docking tag 
      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)

*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE)
c 
      
!     convert LJ atomic parameters to program units  
      do i=1,ntap
        rmin(i)=pnbd1(nbtype(i))/2.d0**(1./6.)
!           abs below is to avoid problems with minus sign in mixing
!           rule
        epsm(i)=abs(pnbd2(nbtype(i))*(1000.d0*4.184/efact)) 
        if(epsm(i).eq.0.d0) rmin(i)=0.1
      end do  


      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      IF(lupdate.and.(.not.linked_cell)) THEN
*--      update shell h neighbor list
        CALL mts_forces(nstep,'u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &       ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &       ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &       ,stressd_p,nnlppf,nnlpp0)
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
         write(kprint,15011) gcpu
         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &        ,ypga,zpga,ngrp,nnlpp0,kprint)
      END IF   
      
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gcpu_u=gcpu
      if(ndocking.eq.0) write(kprint,16011) gcpu
      
*=======================================================================
*---- Zeroes all forces ------------------------------------------------
*=======================================================================
      
      massinfty0 = massinfty*0.5d0
      CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,1,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,prt,utotal,ucns
     &     ,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,fsin14,unb14,cnb14,uslvbon,uslvben,uslvtor,uslvitor
     &     ,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer,enisolv,nstep)
      IF(minimize) THEN
         ok=.FALSE.
         grad_max=-1.0D0
         DO i=1,ntap
            IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
            IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
            IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
         END DO
         
         IF(DABS(grad_max) .LE. eps) THEN
            ok=.TRUE.
            WRITE(kprint,21000) DABS(grad_max)
         ELSE
            if(.not.docking) WRITE(kprint,21500) DABS(grad_max)
         END IF
         
         DO i=1,ntap
           if(mass(i).lt.massinfty0) THEN 
             gpx(i)=-fpx(i)
             gpy(i)=-fpy(i)
             gpz(i)=-fpz(i)
             hpx(i)=gpx(i)
             hpy(i)=gpy(i)
             hpz(i)=gpz(i)
             fpx(i)=hpx(i)
             fpy(i)=hpy(i)
             fpz(i)=hpz(i)
           ELSE
             gpx(i)=0.
             gpy(i)=0.
             gpz(i)=0.
             hpx(i)=0.
             hpy(i)=0.
             hpz(i)=0.
             fpx(i)=0.
             fpy(i)=0.
             fpz(i)=0.
           END IF
         END DO
         
*=======================================================================
*===== MINIMIZATION LOOP ===============================================
*=======================================================================
         
         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         elaps=elapse

         nstep=0
         if(lsasa) THEN 
           lsasa=.false. 
           lsasarun=.true.
         END IF
         DO WHILE(.NOT. ok .AND. nstep .LT. maxstp)
            nstep=nstep+1

            if(nstep.gt.(maxstp-2).and.lsasarun) lsasa=.true. ! last two steps do sasa 
            CALL linmin_total(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &           ,fudgec,xp0,yp0,zp0,xp1,yp1,zp1,fpx,fpy,fpz,fpx1,fpy1
     &           ,fpz1,ntap,fret)
            
            CALL zeroa(fpx,fpy,fpz,ntap,1)
            
c            if(agbnp) call agbnpf(xp0,yp0,zp0,fpx,fpy,fpz,1.d0,enisolv ! call agbnp interface
c     &           ,iret) 
c            write(6,*) " enisolv ", enisolv

            utotal_old=utotal
            CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &           ,1,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,prt,utotal
     &           ,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &           ,coul_bnd_slv_n1,self_slv,fsin14,unb14,cnb14,uslvbon
     &           ,uslvben,uslvtor,uslvitor,uumb,uptors,uitors,ubond
     &           ,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1,coul_bnd_slt
     &           ,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp,eer,enisolv
     &           ,nstep)
            
            CALL add_energies(pme,pressure,slv_exist,slt_exist,ucns,ucos
     &           ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &           ,self_slv,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors
     &           ,uitors,ubond,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1
     &           ,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp,ucosp
     &           ,urcsp,eer,uconf,ucoul,ureal,urecp,pubnd,purecp,puconf
     &           ,pucoul,upconf,upcoul,stressd,stressr,stress_conf
     &           ,stress_tot,co,oc,volume,unitp,press_conf)
            
            grad_max=-1.0D0
            DO i=1,ntap
              if(mass(i).lt.massinfty0) THEN 
                IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
                IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
                IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
              END IF
            END DO

            dutotal=DABS(utotal-utotal_old)
            IF(DABS(grad_max) .LE. eps) THEN 
              if(lsasarun) THEN
                if(lsasa) THEN 
                  ok=.TRUE.
                ELSE
                  lsasa=.true.
                END IF
              ELSE
                ok=.TRUE.
              END IF
            ENDIF
            
            if(docking) THEN 
              if(utotal.lt.utotalmin) utotalmin=utotal
              if(utotal.lt.ucomplex) ucomplex=utotal
            end if
              
            CALL prtacc(pucek,puhyd,puconf,pueng,pucoul,self_slt,fsbond
     &           ,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond,uitors
     &           ,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal
     &           ,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor,uslvitor
     &           ,uumb,rms,temp,temph,tcm,rcm,tempt,tempr,temppr,gr
     &           ,ucepr,stress_tot,press_conf,pressc,press_kin,temppra
     &           ,errca,errhe,errbc,erral,drpca,drpbc,drphe,drpal
     &           ,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_eslvint
     &           ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors,sum_tote
     &           ,sum_ucek,sum_temp,sum_tempt,sum_tempr,sum_temppr
     &           ,sum_temph,sum_rms,sum_pecek,sum_pehyd,sum_peconf
     &           ,sum_pecoul,sum_percip,sum_enb14,sum_ebend,sum_ebond
     &           ,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd,sum_pepot
     &           ,sum_ptote,sum_gr,sum_epcoul,sum_epconf,sum_co,sum_st
     &           ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_volume
     &           ,sum_pv,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd
     &           ,ssm_etotpot,ssm_eslvint,ssm_eebond,ssm_eebend
     &           ,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek,ssm_temp
     &           ,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
     &           ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &           ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &           ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr
     &           ,ssm_epcoul,ssm_epconf,ssm_co,ssm_st,ssm_presst
     &           ,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &           ,ssm_temppra,energy,enisolv,nstep,nstep)
            IF(MOD(nstep,nprint) .EQ. 0) THEN
               WRITE(kprint,20000) DABS(grad_max)
            END IF
            IF(conj_grad .AND. (.NOT. ok)) THEN
               gg=0.0D0
               dgg=0.0D0
               DO i=1,ntap
                 if(mass(i).lt.massinfty0) THEN
                  gg=gg+gpx(i)**2+gpy(i)**2+gpz(i)**2
                  dgg=dgg+(fpx(i)+gpx(i))*fpx(i)+(fpy(i)+gpy(i))
     &                 *fpy(i)+(fpz(i)+gpz(i))*fpz(i)
c$$$               dgg=dgg+fpx(i)**2+fpy(i)**2+fpz(i)**2
                ENDIF
               END DO
               IF(DABS(gg) .EQ. 0.0D0) GOTO 100
               
               gamma=dgg/gg
               
               DO i=1,ntap
                 if(mass(i).lt.massinfty0) THEN 
                   gpx(i)=-fpx(i)
                   gpy(i)=-fpy(i)
                   gpz(i)=-fpz(i)
                   hpx(i)=gpx(i)+gamma*hpx(i)
                   hpy(i)=gpy(i)+gamma*hpy(i)
                   hpz(i)=gpz(i)+gamma*hpz(i)
                   fpx(i)=hpx(i)
                   fpy(i)=hpy(i)
                   fpz(i)=hpz(i)
                 ELSE
                   gpx(i)=0.
                   gpy(i)=0.
                   gpz(i)=0.
                   hpx(i)=0.
                   hpy(i)=0.
                   hpz(i)=0.
                   fpx(i)=0.
                   fpy(i)=0.
                   fpz(i)=0.
                 END IF
               END DO
            END IF
            
*=======================================================================
*---- Write coordinates to a pdb files when required -------------------
*=======================================================================
            
            IF(nascii .NE. 0.and.(.not.docking)) THEN
               IF(MOD(nstep,nascii) .EQ. 0) THEN
                 CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &                ,nprot,protl)
                 CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &                ,zpa)
                 CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &                 ,ypcma,zpcma)
                 CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
                 CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo
     &                ,zpo)
                 fstep=time*DFLOAT(nstep)
                 CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1
     &                ,prsymb,idum)
               END IF
            END IF
            
*=======================================================================
*---- Dump the restart file when required ------------------------------
*=======================================================================
            
            IF(nsave .NE. 0) THEN
               IF(nstep .NE. 0 .AND. MOD(nstep,nsave).EQ.0)THEN
                  grflag= igr .OR. gofr
                  WRITE(kprint,70200)
                  CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co
     &                 ,nprot,protl,mass,tmass,vcax,vcay,vcaz)
                  
                  CALL dumprs(kdump_out,restart_out,nstep,temp,ntap
     &                 ,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,xpg,ypg
     &                 ,zpg,xpcma,ypcma,zpcma,vpx,vpy,vpz,vcax,vcay,vcaz
     &                 ,eta,etap,vh1,3,thermos,sumarray,ssmarray,navg
     &                 ,navg,co,oc,cpress,vco,zz1,lzz,lzz1,grflag,krdf
     &                 ,igrn,nosave)
                  
               END IF
            END IF
         END DO
         
100      CONTINUE

      END IF

*=======================================================================
*---  Print Final Results ----------------------------------------------
*=======================================================================
         
      IF(minimize.and.(.not.docking)) WRITE(kprint,22000)
      CALL add_energies(pme,pressure,slv_exist,slt_exist,ucns,ucos
     &     ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors
     &     ,uitors,ubond,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1
     &     ,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp
     &     ,eer,uconf,ucoul,ureal,urecp,pubnd,purecp,puconf,pucoul
     &     ,upconf,upcoul,stressd,stressr,stress_conf,stress_tot,co
     &     ,oc,volume,unitp,press_conf)
      
      if(.not.docking) nprint=1
      
      CALL prtacc(pucek,puhyd,puconf,pueng,pucoul,self_slt,fsbond
     &     ,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond,uitors
     &     ,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal
     &     ,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor,uslvitor
     &     ,uumb,rms,temp,temph,tcm,rcm,tempt,tempr,temppr,gr
     &     ,ucepr,stress_tot,press_conf,pressc,press_kin,temppra
     &     ,errca,errhe,errbc,erral,drpca,drpbc,drphe,drpal
     &     ,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_eslvint
     &     ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors,sum_tote
     &     ,sum_ucek,sum_temp,sum_tempt,sum_tempr,sum_temppr
     &     ,sum_temph,sum_rms,sum_pecek,sum_pehyd,sum_peconf
     &     ,sum_pecoul,sum_percip,sum_enb14,sum_ebend,sum_ebond
     &     ,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd,sum_pepot
     &     ,sum_ptote,sum_gr,sum_epcoul,sum_epconf,sum_co,sum_st
     &     ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_volume
     &     ,sum_pv,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd
     &     ,ssm_etotpot,ssm_eslvint,ssm_eebond,ssm_eebend
     &     ,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek,ssm_temp
     &     ,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
     &     ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &     ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr
     &     ,ssm_epcoul,ssm_epconf,ssm_co,ssm_st,ssm_presst
     &     ,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &     ,ssm_temppra,energy,enisolv,nstep,nstep)
      
*=======================================================================
*---- Write coordinates to a pdb files and dump restart ----------------
*=======================================================================

      IF(minimize) THEN
         IF(nsave .NE. 0) THEN
            grflag= igr .OR. gofr
            WRITE(kprint,70200)
            CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co
     &           ,nprot,protl,mass,tmass,vcax,vcay,vcaz)
            
            CALL dumprs(kdump_out,restart_out,nstep,temp,ntap,ngrp
     &           ,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,xpg,ypg,zpg,xpcma
     &           ,ypcma,zpcma,vpx,vpy,vpz,vcax,vcay,vcaz,eta,etap,vh1,3
     &           ,thermos,sumarray,ssmarray,navg,navg,co,oc,cpress,vco
     &           ,zz1,lzz,lzz1,grflag,krdf,igrn,nosave)
            
         END IF
         
         IF(nascii .NE. 0.or.docking) THEN
            CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &           ,nprot,protl)
            CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
            CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &           ,ypcma,zpcma)
            CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
            CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo
     &           ,zpo)
            if(docking) THEN 
              fstep=utotalmin
            ELSE
              fstep=time*DFLOAT(nstep)
            END IF
            CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &           ,idum)
         END IF
      END IF 

      if(docking) THEN
        if(idock.eq.0) THEN 
          useparated=utotalmin
        ELSE
          ucomplext=utotalmin
        END IF
        if(idock.gt.0) THEN 
          enerdiss=(useparated-ucomplex)/4.184
          write(kprint,1098) idock,(useparated-ucomplext)/4.184,enerdiss
1098      format(5x,///" ====> Docking round ended ---> ", i5,
     &         " Current dissociation energy(kcal/mol)  =",f10.2,
     &         /5x," ====> Minimal dissociation energy(kcal/mol)  ="
     &         ,f10.2)
        ELSE
          write(kprint,2098) idock,useparated/4.184
2098      format(5x," ====> Docking for separated species round ended -
     &--> ", i5," Final energy(kcal/mol)  =",f10.2)
        END IF
        idock=idock+1
        itry=0
        utotalmin=1.D30
        if(idock.lt.ndocking) THEN 
          write(kprint,11144)
11144     FORMAT(
     &          / '  ===-> Regenerating random ligand coordinates')
!         produce new starting relative ligand-receptor conformation  
1099      call rdocking(ss_index,xp0_s,yp0_s,zp0_s,xpcm,ypcm,zpcm,mass
     &     ,nprot,protl,nstart,nend,jlig,nlig,ntap,dock_radius,xp0,yp0
     &     ,zp0,idock) 
!         check coordinates before stating 
          CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
          call check_coordinates(ntap,co,xpa,ypa,zpa,idum,jdum
     &         ,laux,dummy)
          write(kprint,10789) dummy,idum,jdum
10789     FORMAT(
     &         '  ====> rmin =',f9.4,  ' betweem atoms ', 2i6)
          if(laux.and.itry.lt.1000) THEN 
            itry=itry+1
            write(kprint,*) ' ====> Drug Coordinates are regenerated '
     &           , " -- tentative ", itry
            go to 1099 
          end if
          
!        printouts and restart over 
          CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &                ,nprot,protl)
          CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &         ,zpa)
          CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &         ,ypcma,zpcma)
          CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
          CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo
     &         ,zpo)
          fstep=float(idock)
          CALL plotc(fstep,beta,co,xp0,yp0,zp0,ntap,nres,m1
     &         ,prsymb,idum)
          nprint=1
          CALL prtacc(pucek,puhyd,puconf,pueng,pucoul,self_slt,fsbond
     &         ,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond,uitors
     &         ,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal
     &         ,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor,uslvitor
     &         ,uumb,rms,temp,temph,tcm,rcm,tempt,tempr,temppr,gr
     &         ,ucepr,stress_tot,press_conf,pressc,press_kin,temppra
     &         ,errca,errhe,errbc,erral,drpca,drpbc,drphe,drpal
     &         ,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_eslvint
     &         ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors,sum_tote
     &         ,sum_ucek,sum_temp,sum_tempt,sum_tempr,sum_temppr
     &         ,sum_temph,sum_rms,sum_pecek,sum_pehyd,sum_peconf
     &         ,sum_pecoul,sum_percip,sum_enb14,sum_ebend,sum_ebond
     &         ,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd,sum_pepot
     &         ,sum_ptote,sum_gr,sum_epcoul,sum_epconf,sum_co,sum_st
     &         ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_volume
     &         ,sum_pv,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd
     &         ,ssm_etotpot,ssm_eslvint,ssm_eebond,ssm_eebend
     &         ,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek,ssm_temp
     &         ,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
     &         ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &         ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &         ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr
     &         ,ssm_epcoul,ssm_epconf,ssm_co,ssm_st,ssm_presst
     &         ,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &         ,ssm_temppra,energy,enisolv,nstep,nstep)
          grad_max=-1.0D0
          nprint=100000
          DO i=1,ntap
            if(mass(i).lt.massinfty0) THEN
              IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
              IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
              IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
            ENDIF
          END DO
          WRITE(kprint,20000) DABS(grad_max)
          CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
          goto 1001 ! restart new local minimization
        end if
      end if
      
      if(.not.docking) THEN 
        grad_max=-1.0D0
        DO i=1,ntap
          if(mass(i).lt.massinfty0) THEN
            IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
            IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
            IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
          ENDIF
        END DO
        WRITE(kprint,20000) DABS(grad_max)
        IF(write_grad) CALL prtfrc(kprint,ngrp,grppt,nres,M1,prsymb
     &       ,beta,xp0,yp0,zp0,fpx,fpy,fpz)
      else 
        if(enerdiss.GT.0.d0) THEN
          tds= -3.d0*0.595*log(8.d0*enerdiss/0.595)
        else
          tds=0.d0
        endif
        volume_term = -0.596*log(1661*vol_lig/3.14/vsite**2)
        write(kprint,4098) dock_radius,vsite,enerdiss,tds,volume_term,
     &       enerdiss+tds+volume_term 
4098    FORMAT(//10x,60('='), 
     &       /10x,"=",58x,"=",
     &       /10x,"=",19x,'DOCKING FINAL RESULTS',18x,"="
     &       /10x,"=",58x,"=", 
     &  /10x,"=",'  Eq 32 in  -- > J.Comp.Chem. (2016) 37, 1819-1827',
     &       8x,"=", 
     &       /10x,"=",58x,"=", 
     &       /10x,"=",10x,'Dock radius =',f8.1, ' Vsite =',f9.1,10x,"=",
     &       /10x,"=",58x,"=", 
     &       /10x,"=",10x,'DISSOCIATION ENERGY  =',f10.2,16x,"=",
     &       /10x,"=",10x,'DISSOCIATION ENTROPY =',f10.2,16x,"=",
     &       /10x,"=",10x,'VOLUME TERM          =',f10.2,16x,"=",
     &       /10x,"=",58x,"=", 
     &       /10x,"=",10x,'DISSOCIATION FREE ENERGY = ',f10.2,11x,"=",
     &       /10x,"=",58x,"=", 
     &       /10x,60('=')) 
        
      END IF
*=======================================================================
*---- Compute dynamical matrix -----------------------------------------
*=======================================================================
      
      IF(frequencies.and.(.not.docking)) THEN
         WRITE(kprint,23000) 
         mad=m1
         mbd=n_mat
         CALL comp_dynamic_mat(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &        ,fudgec,xp0,yp0,zp0,mad,mbd,fpx2,fpy2,fpz2,d_mat,wk,eigvl
     &        ,eigvc)
      END IF
      
*=======================================================================
*     Write timing
*=======================================================================
      
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      elaps= -elaps + elapse
      write(kprint,*)
      write(kprint,60030)
      write(kprint,17000) gcpu
      write(kprint,18000) elaps
      IF(nstep .NE. 0) write(kprint,60200) gcpu/DFLOAT(nstep)
      IF(nstep .NE. 0) write(kprint,60300) elaps/DFLOAT(nstep)
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
     &     '      I n t e r m e d i a t e   R e s u l t s        ',
     &     12(' '),'*'/'*',78(' '),'*')
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
13000 FORMAT(/22x,' Temperature has been rescaled ',i5,' times '/)
14000 FORMAT(/22x,'Adjusting bond length to Force Field.'/
     &        22x,'     This will take a while...'//) 
70200 FORMAT('<------ Dumping Restart File ------->'/)
70300 FORMAT(//'<------ Reading Restart File ------->')
70400 FORMAT('<------ Restart File Read in ------->'//)
70100 FORMAT('Velocities have been rescaled ---->'/)

17000 FORMAT(  15x,' Total cpu time for the run     = ',f10.3)
18000 FORMAT(  15x,' Total elapsed time for the run = ',f10.3)
60200 FORMAT(  15x,' Averaged time per step         = ',3x,f7.3)
60300 FORMAT(  15x,' Averaged elapsed per step      = ',3x,f7.3//)
15000 FORMAT(/// 10x,'* * * * r - R E S P A  i s  O N  * * * *'  //)
15500 FORMAT(/   10x,'======= timing ========================='  /)
16011 FORMAT(/5x,'CPUtime for update          =',f10.3)
15011 FORMAT(/5x,'CPUtime for linked cell indexing =',f10.3)
16001 FORMAT(/5x,'CPUtime for h-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
16002 FORMAT(/5x,'CPUtime for l-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
16003 FORMAT(/5x,'CPUtime for m-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
14004 FORMAT(/5x,'CPUtime for n1-contribution  =',f9.4) 
16004 FORMAT(/5x,'CPUtime for n0-contribution  =',f9.4) 
10067 FORMAT(/5x,'THEORIC SPEED UP FOR NON BONDED PART =',f8.2)
10068 FORMAT(/5x,'OVERALL THEORIC SPEED UP =',f8.2/)
20000 FORMAT(
     &     21x,'***************************************'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'/)
21000 FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*        Already at minimum !!        *'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'//)
21500 FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*        Not yet at minimum !!        *'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'//)
22000 FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*        M I N I M I Z A T I O N      *'/
     &     21x,'*                                     *'/
     &     21x,'*       F I N A L   R E S U L T S     *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'//)
23000 FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*     Computing Dynamical Matrix      *'/
     &     21x,'*      and Harmonic Frequencies       *'/
     &     21x,'***************************************'//)
80033 format(/5x,'Expected CPU time for the RUN:',I4,
     &     ' hours and ',I2, ' min',/,  
     &       /5x,' Expected average time per M step:',f8.2,' sec.'/  
     &       /5x,' Expected average time per femto :',f8.2,' sec.'/)  
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
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*     Force Field Parameters for      *'/
     &     21x,'*       the current molecules         *'/
     &     21x,'*                                     *'/
     &     21x,'*                                     *'/
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
1005  FORMAT(80('*'),/ 
     & '*',78(' '),'*',/  
     & '*',30(' '),'DOCKING COMPUTATION',29(' '),'*',/  
     & '*',78(' '),'*',/  
     & '*',10(' '),'LIGAND IS MOLECULE  n. ',i4,' 
     &     WITH',i4,' ATOMS ',20(' '),'*',/    
     & '*',10(' '),'VOLUME OF THE LIGAND IS ',f8.1,' Angs', 31(' '),'*'/,
     & '*',10(' '),'TOTAL NUMBER OF MOLECULES IS ',i4,35(' '),'*'/,
     & '*',10(' '),'DOCKING  RADIUS IS ',f8.2, ' ANGS.',35(' '),'*'/, 
     & '*',10(' '),'DOCKING ROUNDS WILL BE ',i6,39(' '),'*'/, 
     & '*',78(' '),'*',/  
     &     80('*'),/) 
       RETURN
      END

      subroutine rand_rot(natom,x,y,z) 
      implicit NONE 
      integer natom,j
      real*8 x(*),y(*),z(*),axis,pi,cs,ss,theta,duni,xj,yj,zj,dteta
      external duni
      pi = dacos(-1.d0) 
      axis=3*duni()
      dteta=2*pi*duni()
      cs=dcos(dteta)
      ss=dsin(dteta)
      do j=1,natom
c--       asse x
        xj=x(j)
        yj=y(j)
        zj=z(j)
        if (axis.LT.1.0) then
          x(j)=xj
          y(j)=       yj*cs + zj*ss
          z(j)=     - yj*ss + zj*cs
        elseif(axis.GE.1.0.and.axis.LT.2.0) then
          x(j)= xj*cs +       zj*ss
          y(j)=         yj
          z(j)= - xj*ss + zj*cs
        elseif(axis.GE.2.0.and.axis.LT.3.0) then
          x(j)= xj*cs + yj*ss
          y(j)= -xj*ss + yj*cs
          z(j)=                 zj
        endif  
      enddo
      end
      subroutine find_ligand(protl,nprot,atom_prot,jlig)
      implicit none 
      integer protl(*),nprot,m,jlig,mmin,count,j,atom_prot(nprot)
      count=0
      jlig=1
      mmin=10000000
      do  j=1,nprot
        m=protl(1+count)
        atom_prot(j)=m
        count=count+m+1
        if(m.lt.mmin) THEN 
          mmin=m
          jlig=j
        end if
      end do
      
      return
      end
!=========================================================================
      subroutine rdocking(ss_index,xp0_s,yp0_s,zp0_s,xpcm,ypcm,zpcm,mass
     &     ,nprot,protl,nstart,nend,jlig,nlig,ntap,dock_radius,xp0,yp0
     &     ,zp0,idock) 
!=========================================================================
!   produce random relative ligand-receptor configuration for docking 
!=========================================================================

      implicit none 
      REAL*8  xp0_s(*),yp0_s(*),zp0_s(*)
      REAL*8  xp0(*),yp0(*),zp0(*)
      REAL*8  xpcm(*),ypcm(*),zpcm(*),mass(*)
      REAL*8  xpt(ntap),ypt(ntap),zpt(ntap)
      REAL*8  ranf,dummy,dock_radius
      INTEGER ss_index(*),j,jlig,ntap,nstart,nend,protl(*),nprot,nlig
     &     ,idock
      external ranf
      
      dummy=0.0  ! not used 
1099  CALL inicmp(ss_index,xp0_s,yp0_s,zp0_s,xpcm,ypcm,zpcm,mass
     &     ,nprot,protl) 
!         coordinate of the ligand relative to the 0-COM of the ligand

      do j=nstart,nend
        xpt(j)=xp0_s(j)-xpcm(jlig)
        ypt(j)=yp0_s(j)-ypcm(jlig)
        zpt(j)=zp0_s(j)-zpcm(jlig)
      end do
!         three random rotation 
      do j=1,3
        call rand_rot(nlig,xpt(nstart),ypt(nstart)
     &       ,zpt(nstart))
      end do
      do j=1,ntap
        xp0(j)=xp0_s(j)
        yp0(j)=yp0_s(j)
        zp0(j)=zp0_s(j)
      end do
!         now add random COM with docking sphere 
      if(idock.gt.0) THEN 
        xpcm(jlig)=xpcm(jlig)+2.d0*(ranf(dummy)-0.5)*dock_radius
        ypcm(jlig)=ypcm(jlig)+2.d0*(ranf(dummy)-0.5)*dock_radius
        zpcm(jlig)=zpcm(jlig)+2.d0*(ranf(dummy)-0.5)*dock_radius
      ELSE if(idock.eq.0) THEN 
        xpcm(jlig)=xpcm(jlig)+300.d0
        ypcm(jlig)=ypcm(jlig)+300.D0
        zpcm(jlig)=zpcm(jlig)+300.d0
      END IF
      do j=nstart,nend
        xp0(j)=xpt(j)+xpcm(jlig)
        yp0(j)=ypt(j)+ypcm(jlig)
        zp0(j)=zpt(j)+zpcm(jlig)
      end do
      return 
      end
      
      subroutine volume_lig(nstart,nend,rvdw,sasa_nbtype,x,y,z,volume
     &     ,nlig)

      IMPLICIT none
      
      integer nstart,nend, sasa_nbtype(*),nlig
!     args
      real*8 x(*),y(*),z(*),rvdw(*),volume
!     local
      real*8  xl(nlig),yl(nlig),zl(nlig) ! relative to COM ccordinates 
      real*8   xmin,ymin,zmin,xmax,ymax,zmax,rmax,duni,xr,yr,zr,xc,yc,zc
     &     ,vbox,dist
      real*8  rad(nlig)
      integer i,k,npoints,nin
      external duni

      xc=0.d0
      yc=0.d0
      zc=0.d0
      do  i=nstart,nend
        xc = xc + x(i) 
        yc = yc + y(i) 
        zc = zc + z(i) 
      end do
      xc=xc/nlig
      yc=yc/nlig
      zc=zc/nlig
      k=0
      xmax=-1000.0
      ymax=-1000.0
      zmax=-1000.0
      xmin=1000.0
      ymin=1000.0
      zmin=1000.0
      rmax=0.d0
      k=0
      do i=nstart,nend
        k=k+1
        rad(k)=rvdw(sasa_nbtype(i))**2
        xl(k) = x(i)-xc
        yl(k) = y(i)-yc
        zl(k) = z(i)-zc
        if(rvdw(sasa_nbtype(i)).gt.rmax) rmax=rvdw(sasa_nbtype(i)) 
        if(xl(k).gt.xmax) xmax=xl(k)
        if(xl(k).lt.xmin) xmin=xl(k)
        if(yl(k).gt.ymax) ymax=yl(k)
        if(yl(k).lt.ymin) ymin=yl(k)
        if(zl(k).gt.zmax) zmax=zl(k)
        if(zl(k).lt.zmin) zmin=zl(k)
      end do
      xmax=xmax+rmax
      ymax=ymax+rmax
      zmax=zmax+rmax
      xmin=xmin-rmax
      ymin=ymin-rmax
      zmin=zmin-rmax
      vbox=(xmax-xmin)*(ymax-ymin)*(zmax-zmin)
      nin=0
      npoints=vbox/0.2**3  ! grip is 0.2 Angs in  each dir 
      do i=1,npoints
        xr=xmin+duni()*(xmax-xmin)
        yr=ymin+duni()*(ymax-ymin)
        zr=zmin+duni()*(zmax-zmin)
        do k=1,nlig 
          dist=(xl(k)-xr)**2 +(yl(k)-yr)**2 +(zl(k)-zr)**2 ! distance bewteen random point and sphere center
          if(dist.lt.rad(k)) THEN 
            nin=nin+1
            exit
          end if
        end do
      end do
      volume = vbox*dfloat(nin)/dfloat(npoints)
      RETURN
      END 
      
      subroutine gen_itp_file

      use parst
      use cpropar
      use unit 

      IMPLICIT none

*----------------------- LOCAL VARIABLES  --------------------------------
      
      integer i,i1,la,lb,lc,ld,im,indg(m1)
      real*8  phase,sign,gcpu
      gcpu = unite*avogad/1000.d0   ! conversion orac units -> kJ mol-1

*--------------------EXECUTABLE STATEMENTS  ----------------------------------


c============ [defaults ] directive ===================================

      write(kitp,'(A,/)') 
     &  "; Set combination rules for NB parameters and fudge factors
     & for 14 interactions"
      write(kitp,'(A)') 
     &     "[ defaults ]"
      write(kitp,'(A)') 
     &     "; NBFUNC    COMB-RULE       GEN-PAIRS  FUDGELJ      FUDGEQQ"
      write(kitp,'(A)') 
     &     "1           2               YES        0.5000       0.8333"

c============ [ atomstypes ] directive ===================================
      write(kitp,'(//,A)') 
     &  "; Atomic types, masses and charges (au) are defined below"
      write(kitp,'(A)') 
     &  "; N.B.: use maxwarn 100 to ignore 'Overriding atomtypes'  "
      write(kitp,'(A,/)') 
     &  "; N.B: SIGMA is in nm and EPSILON in kj mol-1             "

      write(kitp,'(A)') 
     &     "[ atomtypes ]"
      write(kitp,'(A)')
     &     ";NAME   AT.NUM  MASS     CHARGE  PTYPE  SIGMA   EPSILON"

      im=6   ! defauls is all carbon 
      do i=1,ntap
        if(betb(i)(1:1).eq."h") im=1
        if(betb(i)(1:1).eq."c") im=6
        if(betb(i)(1:1).eq."o") im=8
        if(betb(i)(1:1).eq."n") im=7
        if(betb(i)(1:1).eq."p") im=15
        if(betb(i)(1:1).eq."s") im=16
        if(betb(i)(1:1).eq."f") im=9
        if(betb(i)(1:1).eq."i") im=53
        if(betb(i)(1:2).eq."cl") im=17
        if(betb(i)(1:2).eq."br") im=35
        write(kitp,'(A6,i6,2f10.3,"   A  ",2f9.4 )') betb(i),im,mass(i)
     &       ,chrge(i)*dsqrt(unitc),0.2*pnbd1(nbtype(i))/(2**(1./6.)) 
     &       ,pnbd2(nbtype(i))*4.184
      end do

      write(kitp,'(//,A)') 
     &  "; Molecule topology/parameters starts below"

      write(kitp,'(A)') 
     &     "[ moleculetype ]"
      write(kitp,'(A)') 
     &     "; Name               nrexcl"
      write(kitp,'(A)') 
     &     "  name-p                  3 "

c============ [ atoms ] directive ===================================

      write(kitp,'(//,A,/)') 
     &  "; Atomic types, pdb names  and groups are defined below" 
      write(kitp,'(A)') 
     &     "[ atoms ]"
      do i=1,ngrp
        do i1=grppt(1,i),grppt(2,i)
          indg(i1)=i    ! group index of atom i1
        end do
      end do

      write(kitp,'(A)')
     &";AT.NUM   TYPE   RESID   RESNAME PDB-NAME  IGRP    CHRGE"
     &,"    MASS" 
      do i=1,ntap
        write(kitp,'(I6,5x,A6,3x,"1",6x,"LIG",5x,A6,I4,2f11.4)') i
     &       ,betb(i),beta(i),indg(i),chrge(i)*dsqrt(unitc),mass(i)
      end do

c============ [ bonds  ] directive ===================================

      write(kitp,'(//,A,/)') 
     &  "; Stretching parameters are defined below" 
      write(kitp,'(A)') 
     &     "[ bonds ]"
      write(kitp,'(A)')
     &     ";   I1    I2  TYPE      r0(nm)   K(kJ/nm^2/mol) " 
      do i=1,lbond
        la=lbnd(1,i)
        lb=lbnd(2,i)
        write(kitp,'(2I6,5x,"1", f12.4,f15.3)') la,lb,potbo(i,2)*0.1
     &       ,potbo(i,1)*gcpu*100.d0
      end do

c============ [ angles   ] directive ===================================

      write(kitp,'(//,A,/)') 
     &  "; Bending parameters are defined below" 
      write(kitp,'(A)') 
     &     "[ angles ]"
      write(kitp,'(A)')
     &     ";   I1    I2    I3  TYPE    theta(deg)  K(kJ/rad^2/mol)" 

      do i=1,lbend
        la=lbndg(1,i)
        lb=lbndg(2,i)
        lc=lbndg(3,i)
        write(kitp,'(3I6,5x,"1", f12.2,f15.5)') la,lb,lc,potbe(i,2)
     &       *180.d0/pi,potbe(i,1)*gcpu
      end do

c============ [ dihedrals ] directive (proper) ==============================

      write(kitp,'(//,A,/)') 
     &  "; Proper Torsion parameters are defined below" 
      write(kitp,'(A)') 
     &     "[ dihedrals ]"
      write(kitp,'(A)')
     &   ";   I1    I2    I3    I4"
     &   ,"  TYPE   theta(deg) K(kJ/rad^2/mol)     N" 

      do i=1,ltors 
        la=ltor(1,i)
        lb=ltor(2,i)
        lc=ltor(3,i)
        ld=ltor(4,i)
        sign=potto(i,2)/abs(potto(i,2)) ! handle negative V0
        if(potto(i,1).gt.0.d0) THEN ! sign of v0 is phase 180/0
          phase=180.0 
        ELSE
          phase=0.0 
        END IF
        write(kitp,'(4I6,5x,"1", 2f12.4,i12)') la,lb,lc,ld,phase, sign
     &       *abs(potto(i,1))*gcpu,int(abs(potto(i,2)))
      end do


c============ [ dihedrals ] directive (improper) ==============================

      write(kitp,'(//,A,/)') 
     &  "; ImProper Torsion parameters are defined below" 
      write(kitp,'(A)') 
     &     "[ dihedrals ]"
      write(kitp,'(A)')
     &   ";   I1    I2    I3    I4"
     &   ,"  TYPE   theta(deg) K(kJ/rad^2/mol)     N" 

      do i=1,litor 
        la=litr(1,i)
        lb=litr(2,i)
        lc=litr(3,i)
        ld=litr(4,i)
        sign=potit(i,2)/abs(potit(i,2)) ! handle negative V0
        if(potto(i,1).gt.0.d0) THEN ! sign of v0 is phase 180/0
          phase=180.0 
        ELSE
          phase=0.0 
        END IF
        write(kitp,'(4I6,5x,"1", 2f12.4,i12)') la,lb,lc,ld,phase, sign
     &       *abs(potit(i,1))*gcpu,int(abs(potit(i,2)))
      end do
      return 
      end 
