************************************************************************
*   Time-stamp: <98/03/20 16:22:12 marchi>                             *
*                                                                      *
*   Compute averages and do some analysis at time step M               *
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 31 1998 -                                     *
*                                                                      *
************************************************************************


*=======================================================================
*---        Computes all K.E.s, temperatures 
*=======================================================================
            
      CALL kinetic(ss_index,co,nato_slv,nmol,cnstpp_slv,ntap
     &     ,cnstpp,vpx,vpy,vpz,vcax,vcay,vcaz,nprot,protl,mass
     &     ,tmass,wmtp,ucek,pucek,temp,tempt,tempr,temppr,tcm
     &     ,rcm,stresstk,cpress,isostress,vco,masspp,ucepr
     &     ,temppra,massinfty)
      
*=======================================================================
*---        Gather all energy terms 
*=======================================================================
            
            
*--         Energy terms for the solvent         -----------------------
            
1010  CONTINUE
      fsrtalc=0.0D0
      IF(slv_exist) THEN
        uconf=ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m+urcs_h +urcs_l
     &       +urcs_m+coul_bnd_slv+conf_bnd_slv_n1+coul_bnd_slv_n1
     &       +self_slv+uumb
        ucoul=ucos_h+ucos_l+ucos_m+urcs_h+urcs_l+urcs_m+coul_bnd_slv
     &       +coul_bnd_slv_n1+self_slv
        ureal=ucos_h+ucos_l+ucos_m+coul_bnd_slv_n1
        urecp=urcs_h+urcs_l+urcs_m+self_slv+coul_bnd_slv
      END IF
            
*---        Energy terms for the protein         -----------------------
            
      IF(slt_exist) THEN 
        pubnd=uptors+uitors+ubond+ubend
        unbond=ucnp_m + ucnp_l + ucnp_h
        cnbond=ucop_m + ucop_l + ucop_h
        purecp=urcp_h+urcp_l+urcp_m
        puconf= unbond + conf_bnd_slt_n1 + uumb
        pucoul= cnbond + purecp + coul_bnd_slt + coul_bnd_slt_n1
     &       + self_slt
      END IF
            
            
*---        Mixed terms                          -----------------------
            
      IF(slv_exist .AND. slt_exist) THEN 
        upconf= ucnsp_h + ucnsp_l + ucnsp_m + 
     &       ucosp_h + ucosp_l + ucosp_m + 
     &       urcsp_h + urcsp_l + urcsp_m
        upcoul= ucosp_h + ucosp_l + ucosp_m + 
     &       urcsp_h + urcsp_l + urcsp_m
      END IF
      if(addstrcom) pubnd=pubnd+ucom
            
      IF(pressure) THEN
        DO i=1,3
          DO j=1,3
            stressd(i,j)=stressd_h(i,j)+stressd_l(i,j)+stressd_m(i,j)
     &         +  stresst0(i,j)
            stressr(i,j)=stressr_h(i,j)+stressr_l(i,j)+stressr_m(i,j)
          END DO
        END DO

        CALL comp_stress_conf(stressd,stressr,stress_conf,oc,volume,
     &       unitp,press_conf)
        CALL comp_stress_kinetic(vcax,vcay,vcaz,tmass,co,nprot,volume
     &       ,unitp,stress_kin,press_kin)
        CALL dcopy(9,stress_kin,1,stress_tot,1)
        CALL daxpy(9,1.0D0,stress_conf,1,stress_tot,1)
      END IF
            
*----       If pme=T no solute-solvent separation ----------------------
            
      IF(pme) THEN 
        IF(slv_exist.AND.(.NOT.slt_exist)) THEN
          uconf=uconf + eer_h + eer_l + eer_m
          ucoul=ucoul + eer_h + eer_l + eer_m
          urecp= eer_h + eer_l + eer_m
        END IF
        IF(slt_exist.AND.(.NOT.slv_exist)) THEN
          purecp=eer_h + eer_l + eer_m
          pucoul=pucoul+eer_h + eer_l + eer_m
        END IF
        IF(slt_exist.AND.slv_exist) THEN 
          upconf=upconf + eer_h + eer_l + eer_m
          upcoul=upcoul + eer_h + eer_l + eer_m
        END IF
      END IF

!#ifdef DEBUGALCHEMY
c      write(80,8097)  "Slv", rtime/1000.d0,ucoul,ucos_h,ucos_l,ucos_m
c     &     ,coul_bnd_slv,coul_bnd_slv_n1,self_slv
c      write(80,8097)  "Slt", rtime/1000.d0,pucoul*1.d4
c     &     ,ucop_m*1.d4,coul_bnd_slt*1.d4,coul_bnd_slt_n1*1.d4,self_slt
c     &     *1.D4
c      write(81,8099) "Slt", falch_slt*1.d4,coul_bnd_slt_n1*1.d4
c     &     ,cnbond*1.d4,dwrk_al*1.d4
c8099  format(A5,2x,4f12.5)  
c      write(80,8097)  "Mix", rtime/1000.d0,upcoul*1.d4,ucosp_h*1.d4,
c     &     ucosp_l*1.d4,ucosp_m*1.d4,eer_l*1.d4
c8097  format(A3,1x,f12.4,10G15.5)
!#endif
      IF(energy_then_die) THEN 
!       print configuration and corresponding energy then exit
        IF(ascii_nocell) THEN
          CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xpo,ypo,zpo)
        ELSE
          CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
          CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
        END IF
        fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
        para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
        CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &       ,para0)
        go to 1020
      END IF
*=======================================================================
*---   Now for each timestep of shell M do inline analysis, scale    ---
*---     temperatures, compute averages and dump hystory files       ---
*=======================================================================
            
            
*=======================================================================
*---- Compute instantaneous X-rms --------------------------------------
*=======================================================================
            
      IF(anxrms) CALL calc_xrms(anxca,anxbc,anxhe,anxal,anprot,annpro
     &     ,anpoint,protl,wca,whe,wbc,xp0,yp0,zp0,xpt0,ypt0,zpt0
     &     ,nato_slt,errca,errhe,errbc,erral,drpca,drpbc,drphe,drpal)
            
*=======================================================================
*---- Compute averaged structure ---------------------------------------
*=======================================================================
            
      IF(avg_str) THEN
        IF(avg_ca) CALL calc_avg_str(anxrms_cell,protl,wca,xpt0,ypt0
     &       ,zpt0,xp_avg,yp_avg,zp_avg,qt,xp0,yp0,zp0,nato_slt
     &       ,iter_avg)
        IF(avg_he) CALL calc_avg_str(anxrms_cell,protl,whe,xpt0,ypt0
     &       ,zpt0,xp_avg,yp_avg,zp_avg,qt,xp0,yp0,zp0 ,nato_slt
     &       ,iter_avg)
      END IF
            
      IF(navg_str .NE. 0) THEN
        IF(MOD(ninner,navg_str) .EQ.0) THEN
          CALL dcopy(nato_slt,xp_avg,1,xpo,1)
          CALL dcopy(nato_slt,yp_avg,1,ypo,1)
          CALL dcopy(nato_slt,zp_avg,1,zpo,1)
          fact=1.0D0/DFLOAT(iter_avg)
          CALL dscal(nato_slt,fact,xpo,1)
          CALL dscal(nato_slt,fact,ypo,1)
          CALL dscal(nato_slt,fact,zpo,1)
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          IF(avg_ca) WRITE(kavg,80000) 
          IF(avg_he) WRITE(kavg,80100) 
          CALL plotd(fstep,kavg,beta,xpt0,ypt0,zpt0,xpo,ypo,zpo
     &         ,nato_slt,nres,m1,prsymb)
        END IF
      END IF
      IF(navg_str_xrms .NE. 0) THEN
        IF(MOD(ninner,navg_str_xrms) .EQ.0) THEN
          CALL dcopy(nato_slt,xp_avg,1,xpo,1)
          CALL dcopy(nato_slt,yp_avg,1,ypo,1)
          CALL dcopy(nato_slt,zp_avg,1,zpo,1)
          fact=1.0D0/DFLOAT(iter_avg)
          CALL dscal(nato_slt,fact,xpo,1)
          CALL dscal(nato_slt,fact,ypo,1)
          CALL dscal(nato_slt,fact,zpo,1)
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          CALL calc_avg_xrms(avg_ca,avg_he,fstep,kavg_xrms,xpt0,ypt0
     &         ,zpt0,xpo,ypo,zpo,wca,whe,wbc,protl,nato_slt)
        END IF
      END IF
            
            
*=========== Write instantaneous results, averages and fluctuations ====
            

1020  CALL prtacc(pucek,puhyd,puconf,pueng,pucoul
     &     ,self_slt,fsbond,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond
     &     ,uitors,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal
     &     ,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor,uslvitor,uumb,rms
     &     ,temp,temph,tcm,rcm,tempt,tempr,temppr,gr,ucepr,stress_tot
     &     ,press_conf,pressc,press_kin,temppra,errca,errhe,errbc,erral
     &     ,drpca,drpbc,drphe,drpal,sum_econf,sum_ecoul,sum_enbnd
     &     ,sum_etotpot,sum_eslvint,sum_eebond,sum_eebend,sum_eeptors
     &     ,sum_eeitors,sum_tote,sum_ucek,sum_temp,sum_tempt,sum_tempr
     &     ,sum_temppr,sum_temph,sum_rms,sum_pecek,sum_pehyd,sum_peconf
     &     ,sum_pecoul,sum_percip,sum_enb14,sum_ebend,sum_ebond
     &     ,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd,sum_pepot,sum_ptote
     &     ,sum_gr,sum_epcoul,sum_epconf,sum_co,sum_st,sum_presst
     &     ,sum_press,sum_pressc,sum_pressk,sum_volume,sum_pv
     &     ,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot
     &     ,ssm_eslvint,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors
     &     ,ssm_tote,ssm_ucek,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr
     &     ,ssm_temph,ssm_rms,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul
     &     ,ssm_percip,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &     ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul
     &     ,ssm_epconf,ssm_co,ssm_st,ssm_presst,ssm_press,ssm_pressc
     &     ,ssm_pressk,ssm_volume,ssm_pv,ssm_temppra,energy,enisolv
     &     ,ninner,nstep)
      IF(energy_then_die) then
#ifdef _MPI_      
      CALL MPI_Barrier(MPI_COMM_WORLD,iret) 
      write(kprint,10034)
      call MPI_FINALIZE(ierr)
#else        
      write(kprint,10034)
#endif
10034 FORMAT(
     &     10x,"======================================================"/
     &     10x,"=                                                    ="/ 
     &     10x,"=            STOP:    energy_then die                ="/ 
     &     10x,"=                                                    ="/ 
     &     10x,"======================================================")  
        STOP
      endif
*=========== Dump files ===============================================
      IF(ndipole.gt.0) THEN
        IF(MOD(ninner,ndipole).EQ.0)THEN
          CALL comp_dip2(xp0,yp0,zp0,chrge,dips,ntap)
          dip_tot=dip_tot+dips(1)**2+dips(2)**2+dips(3)**2
          ndip=ndip+1 
          aux = dip_tot*4.d0*pi*efact/(3.d0*volume*gascon*temp
     &         *float(ndip))
          do j=1,3
            dips(j) = 4.8*dips(j)*dsqrt(unitc)
          end do
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          write(kdipole,106) fstep, (dips(j),j=1,3),aux
106       FORMAT(f12.3, 3e15.5,f10.3)
        END IF
      END IF
#ifdef _MPI_
      if(iproc.eq.0) then
        if(mod(ninner,nprint).eq.0) then
          write(*,1031) mstep * time
1031      FORMAT( 15x,'==== current simulation time = ', f12.1,2i10)  
        endif
      endif
#endif
      IF(nascii .NE. 0) THEN
        IF(MOD(ninner,nascii) .EQ. 0) THEN
          IF(ascii_nocell) THEN
            CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xpo,ypo,zpo)
          ELSE
            CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
            CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
          END IF
!         call print_cm(xp0,yp0,zp0,mass,ntap
!    &         ,massinfty)
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
          CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &         ,para0)
        END IF
      END IF
      if(time_steer.and.alchemy) THEN 
        if(nplotalch0.gt.0) THEN 
          if(MOD(ninner,nplotalch0).EQ.0) THEN 
            if(( (.not.lrm_end) .or. (.not.lrmq_end) ) .or.
     &           ( (.not.ladd_end) .or. (.not.laddq_end) ) )THEN
              call print_alchemic_wrk(rtime,ladd,laddq,lrm,lrmq
     &             ,wrk_alchemy,kplot_alch0,efact)
            end if
          endif
        END IF
      END IF
*=========== Dump configurations to a .PDB file =========================

      IF(ndcd /= 0) THEN
        IF(MOD(ninner,ndcd) == 0)THEN
          CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
          CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
          call print_dcd(ntap,co,xpo,ypo,zpo,beta)
        END IF
      END IF
            
      IF(nplot.GT.0) THEN
        IF(MOD(ninner,nplot).EQ.0)THEN
          CALL mts_plotp(beta,mback,nbone,xp0,yp0,zp0,ntap,nres,m1
     &         ,prsymb)
        END IF
      END IF
          
      IF( nplot_fragm .GT. 0 ) THEN
        IF( MOD(ninner,nplot_fragm).EQ.0 ) THEN
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
! print a file containing num. of atoms time, the co matrix and
! coordinates
          if ( print_unformatted ) THEN
            write (kplot_fragm) ntot_fragm
#ifdef _MPI_
            if ( rem_run .or. sim_tempering ) then
#else
            if ( sim_tempering ) then
#endif
              write (kplot_fragm) para_index, sngl(fstep),
     &             sngl(2.d0*co(1,1)), sngl(2.d0*co(1,2)),
     &             sngl(2.d0*co(1,3)), sngl(2.d0*co(2,1)),
     &             sngl(2.d0*co(2,2)), sngl(2.d0*co(2,3)),
     &             sngl(2.d0*co(3,1)), sngl(2.d0*co(3,2)),
     &             sngl(2.d0*co(3,3))
            else
              write (kplot_fragm) sngl(fstep), sngl(2.d0*co(1,1)),
     &             sngl(2.d0*co(1,2)), sngl(2.d0*co(1,3)),
     &             sngl(2.d0*co(2,1)), sngl(2.d0*co(2,2)),
     &             sngl(2.d0*co(2,3)), sngl(2.d0*co(3,1)),
     &             sngl(2.d0*co(3,2)), sngl(2.d0*co(3,3))
            endif
          else
            write (kplot_fragm,*) ntot_fragm
#ifdef _MPI_
            if ( rem_run .or. sim_tempering ) then
#else
            if ( sim_tempering ) then
#endif
              para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
              write (kplot_fragm,11) fstep,
     &             2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &             2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &             2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3),
     &             para_index,para0

11            format('Time:',1x,f15.2,1x,'CO matrix:',1x,9f10.5,1x,
     &             'Ensemble:',1x,2i4)
            else
              write (kplot_fragm,10) fstep,
     &             2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &             2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &             2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3)
10            format('Time:',1x,f15.2,1x,'CO matrix:',1x,9f10.5)
            endif
          endif
          if(.not.xyz_nocell) THEN 
            CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xpo,ypo,zpo)
          ELSE
            CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
            CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
          ENDIF
            do i = 1, nfragm
            fragm_1 = fragm(1,i)
            fragm_2 = fragm(2,i)
            CALL mts_plot_fragm(fragm_1,fragm_2,beta,xpo,ypo,zpo
     &           ,ntap)
          end do
        ENDIF
      END IF
            
      IF(fragm_dist) THEN
        IF(MOD(ninner,nfragm_dist) .EQ. 0) THEN
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          CALL write_fragm_dist(fstep,kfragm_dist,fragm
     &         ,nfragm,xpa,ypa,zpa,co)
        END IF
      END IF

      IF(nplot_center .NE. 0) THEN
        IF(MOD(ninner,nplot_center) .EQ. 0) THEN
          CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
          CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
! print a file containing num. fo atoms time, the co matrix and coordinates
          IF ( print_unformatted ) THEN
            write (kplot_center) ntap
#ifdef _MPI_
            if ( rem_run .or. sim_tempering ) then
#else
            if ( sim_tempering ) then
#endif
              write (kplot_center) para_index, sngl(fstep),
     &             sngl(2.d0*co(1,1)), sngl(2.d0*co(1,2)),
     &             sngl(2.d0*co(1,3)), sngl(2.d0*co(2,1)),
     &             sngl(2.d0*co(2,2)), sngl(2.d0*co(2,3)),
     &             sngl(2.d0*co(3,1)), sngl(2.d0*co(3,2)),
     &             sngl(2.d0*co(3,3))
            else
              write (kplot_center) sngl(fstep), sngl(2.d0*co(1,1)),
     &             sngl(2.d0*co(1,2)), sngl(2.d0*co(1,3)),
     &             sngl(2.d0*co(2,1)), sngl(2.d0*co(2,2)),
     &             sngl(2.d0*co(2,3)), sngl(2.d0*co(3,1)),
     &             sngl(2.d0*co(3,2)), sngl(2.d0*co(3,3))
            endif
          ELSE
            write (kplot_center,*) ntap
#ifdef _MPI_
            if ( rem_run .or. sim_tempering ) then
#else
            if ( sim_tempering ) then
#endif
              write (kplot_center,11) para_index, fstep,
     &             2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &             2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &             2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3)
            else
              write (kplot_center,10) fstep,
     &             2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &             2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &             2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3)
            endif
          ENDIF
          CALL plot_center(fstep,beta,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &           ,chrge)
        END IF
      END IF
            
*=========== Dump the hystory of the run when required =================
            
      IF(nconf.GT.0) THEN
        IF(MOD(ninner,nconf).EQ.0 .AND. ninner.GT.mrject)THEN
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          CALL write_confc(co,xp0,yp0,zp0,ntap,fstep
     &         ,ninner,nconf,divide_records,atom_record)
        END IF
      END IF
            
            
      IF(gofr) THEN
        IF(MOD(ninner,gofr_ncomp).EQ.0 .AND. nstep.GT.nrject) THEN
          l2=maxint
          CALL calc_gofr(nato_slt,nato_slv,type_slv,ss_index,atomp
     &         ,gofr_neighbor,gofr_intra,co,xpa,ypa,zpa,xpga,ypga,zpga
     &         ,wca,whe,delrg,ntap,ngrp,grppt,l2,krdf,ngrdon,gofr_cut
     &         ,nnlpp0)
          IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
        END IF
      END IF
            
      IF(gofr) THEN
        IF(MOD(ninner,gofr_nprint).EQ.0 .AND. nstep.GT.nrject)THEN
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          vol_gofr=sum_volume/DFLOAT(ninner)
          IF(slt_exist) THEN
            offset=0
            CALL write_gofrp(.NOT.gofr_avg,fstep,krdf
     &           ,maxint,offset,wca,whe,nato_slt,delrg,gofr_cut,ngrdon
     &           ,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
          END IF
          IF(slv_exist .AND. slt_exist) THEN
            CALL write_gofrw(.FALSE.,fstep,krdf,maxint
     &           ,g1,delrg,gofr_cut,itype_slv,betab_slv,vol_gofr
     &           ,ntype_slv,nmol,ngrdon,3,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
          ELSE IF(slv_exist .AND. (.NOT. slt_exist)) THEN
            CALL write_gofrw(.TRUE.,fstep,krdf,maxint,g1
     &           ,delrg,gofr_cut,itype_slv,betab_slv,vol_gofr,ntype_slv
     &           ,nmol,ngrdon,3,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
          END IF
        END IF
      END IF
            
*=======================================================================
*---- Initialize G of Rs when required ---------------------------------
*=======================================================================
            
      IF(gofr .AND. gofr_avg  .AND. nstep .GT. nrject) THEN
        IF(MOD(ninner,gofr_navg) .EQ. 0) THEN
          CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
        END IF
      END IF
            
*=======================================================================
*----- Write topology --------------------------------------------------
*=======================================================================
            
      IF(prttopl) THEN
        IF(MOD(ninner,ntop_print) .EQ. 0) THEN
          fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
          IF(top_bonds(1) .GT. 0) CALL
     &         write_bonds(ktopol,fstep,top_bonds,lbnd,lbond,xp0,yp0
     &         ,zp0)
          IF(top_bendings(1) .GT. 0) CALL
     &         write_bends(ktopol,fstep,top_bendings,lbndg,lbend,xp0,yp0
     &         ,zp0)
          IF(top_ptors(1) .GT. 0) CALL write_tors('P'
     &         ,ktopol,fstep,top_ptors,ltor,ltors,xp0,yp0,zp0)
          IF(top_itors(1) .GT. 0) CALL write_tors('I'
     &         ,ktopol,fstep,top_itors,litr,litor,xp0,yp0,zp0)
        END IF
      END IF
            
*=========== Scale the temperature when required =======================
            
      IF(ninner-mrject.EQ.0)THEN
               
        CALL kinetic(ss_index,co,nato_slv,nmol,cnstpp_slv,ntap,cnstpp
     &       ,vpx,vpy,vpz,vcax,vcay,vcaz,nprot,protl,mass,tmass,wmtp
     &       ,ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk
     &       ,cpress,isostress,vco,masspp,ucepr,temppra,massinfty)
        
        IF(slv_exist) ustot=(ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m
     &       +urcs_h+urcs_l+urcs_m+uslvbon+uslvben+uslvtor+uslvitor
     &       +coul_bnd_slv+self_slv)*efact/1000.0D0
               
        uptot=0.d0
        IF(slt_exist) uptot=(urcp_h+urcp_l+urcp_m+ucop_h+ucop_l+ucop_m
     &       +ucnp_h+ucnp_l+ucnp_m+ubond+ubend+uptors+uitors
     &       +coul_bnd_slt+self_slt)*efact/1000.d0
        ektot = (ucek+pucek)/1000.d0
        nscal=0
               
        upstot=0.d0
        IF(slv_exist .AND. slt_exist) upstot = (urcsp_h+urcsp_l+urcsp_m
     &       +ucosp_h+ucosp_l +ucosp_m+ucnsp_h+ucnsp_l+ucnsp_m+eer_m
     &       +eer_l+eer_h)*efact/1000.d0
        E0 = ektot+ustot+uptot+upstot
        nrject=-1
        mrject=-1
        ninner = 0
        nstep=0
        write(kprint,13000) nscale
      END IF
      aux=efact/1000.d0  
!     write(6,1083) rtime,aux*coul_bnd_slt,aux*ucop_m,aux*cnb14,aux
!     &     *self_slt
!1083 format(f10.1," Erf", f10.3, " Direct", 2f10.3, " Self", f10.3)
      IF(annealing .AND. ninner-mrject .GT. 0) THEN
        CALL anneal(annealing_fact,vpx,vpy,vpz,ntap)
        IF(cnstpp .NE. 0) THEN
          CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap,cnstp,dssp
     &         ,coeffp,cnstpp,mass,dnit,cnst_protp,cnst_protl,mim_lim
     &         ,iret,errmsg)
          IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
        END IF
        CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &       ,vcaz)
      END IF
            
      IF(ninner-mrject .LT. 0) THEN
              
        IF(DABS(t-temp).GT.dtemp .OR. MOD(ninner,scale).EQ.0) THEN
          nscale=nscale+1
          CALL change_origin(-1,nprot,protl,vpx,vpy,vpz,vpx,vpy,vpz
     &         ,vcax,vcay,vcaz,co)
          CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,.TRUE.
     &         ,massinfty)
          IF(cnstpp .NE. 0) THEN
            CALL rattle_correc(tm,xp0,yp0,zp0,vpx,vpy,vpz,ntap,cnstp
     &           ,dssp,coeffp,cnstpp,mass,dnit,cnst_protp,cnst_protl
     &           ,mim_lim,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
          END IF
          CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &         ,vcaz)
          IF(.not.start_conf) write(kprint,70100)
        END IF
        IF(cpress) THEN
          IF(DABS(t-temppra).GT.dtemppr .OR. MOD(ninner,scale).EQ.0)THEN
            CALL set_tempp(masspr,vco,temppra,t)
            IF(isostress) THEN
              CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,nboxcnst
     &             ,iret,errmsg)
              IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF
                  
            IF(.not.start_conf) write(kprint,70120)
          END IF
        END IF
        IF(thermos) THEN
          IF(DABS(t-temph).GT.dtemph .OR. MOD(ninner,scale).EQ.0) THEN
            CALL set_tempt(neta,qmass,etap,temph,t)
            write(kprint,70130)
          END IF
        END IF
               
        ig=0
      ENDIF
      IF(always_scale_nose) THEN
        IF(DABS(t-temph) .GT. dtemph .OR. MOD(ninner,scale) .EQ. 0) THEN
          CALL set_tempt(neta,qmass,etap,temph,t)
          write(kprint,70130)
        END IF
        IF(cpress) THEN
          IF(DABS(t-temppra).GT.dtemppr .OR. MOD(ninner,scale).EQ.0)THEN
            CALL set_tempp(masspr,vco,temppra,t)
            IF(isostress) THEN
              CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,nboxcnst
     &             ,iret,errmsg)
              IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF
            write(kprint,70120)
          END IF
        END IF
      END IF
c--   frequency-based printouts of target replica in cycling SGE
      if(always_accept) THEN 
        IF(para_index.eq.target_replica) THEN
          n_always_accept=n_always_accept+1
          if(ninner.GT.mrject) THEN  
          if(mod(n_always_accept,nfreqsge).eq.0) THEN
            if(print_always.eq.1) THEN 
              CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
              CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
              fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
              para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
              CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1,prsymb
     &             ,para0)
            ELSE IF (print_always.eq.2) THEN 
              if(ntot_fragm.eq.0) THEN
                errmsg=
     &        "No fragment defined: coordinates will not be printed  "
                CALL xerror(errmsg,80,1,2)
                STOP
              endif
              write (kplot_fragm,*) ntot_fragm
              fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
              para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
              write (kplot_fragm,11) fstep,
     &             2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &             2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &             2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3),
     &             para_index,para0
              
              CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl)
              CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo,zpo)
              do i = 1, nfragm
                fragm_1 = fragm(1,i)
                fragm_2 = fragm(2,i)
                CALL mts_plot_fragm(fragm_1,fragm_2,beta,xpo,ypo,zpo
     &               ,ntap)
              end do
            ENDIF
          endif
          endif
        END IF    
      end if
c--   check whether alchemical trsfs are over 
      if(.not.alchemy) go to 1009 
      if(near0(ladd).and.(.not.ladd_end).and.nat_added.ne.0)
     &     THEN
        write(kprint,10985) (1.d0-ladd)*100.d0,rtime
10985   FORMAT("!! LJ-Alchemical transformation (add)  ended:", 
     &       f8.2,"% done at total simulation time",f10.2  )
        ladd_end=.true. 
      ENDIF
      if(near0(laddq).and.(.not.laddq_end).and.nat_added.ne.0)
     &     THEN
        write(kprint,10993) (1.d0-laddq)*100.d0,rtime
10993   FORMAT("!! QQ-Alchemical transformation (add)  ended:", 
     &       f8.2,"% done at total simulation time",f10.2  )
        laddq_end=.true. 
      ENDIF
      aux=lrm+lambda0(atom_removed(1))
      if(aux.lt.0.and.(.not.lrm_end).and.nat_removed.ne.0)THEN
        write(kprint,10984)  lrm*100.d0,rtime
10984   FORMAT("!! LJ-Alchemical transformation (remove) ended:", 
     &       f8.2,"% done at total simulation time",f10.2  )
        lrm_end=.true. 
!               compute the work due to recp of ghost atoms (see
!               initialize_wrk) 
      END IF
      aux=lrmq+lambdaq0(atom_removed(1))
      if(aux.lt.0.and.(.not.lrmq_end).and.nat_removed.ne.0)THEN
        write(kprint,10994)  lrmq*100.d0,rtime
10994   FORMAT("!! QQ-Alchemical transformation (remove) ended:", 
     &       f8.2,"% done at total simulation time",f10.2  )
        lrmq_end=.true. 
!               compute the work due to recp of ghost atoms (see
!               initialize_wrk) 
      ENDIF
1009  CONTINUE



