***********************************************************************
*   Time-stamp: <98/03/20 19:56:59 marchi>                             *
*                                                                      *
*   Dump restart file and do tests at time H                           *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 31 1998 -                                     *
*                                                                      *
************************************************************************
 
*=======================================================================
*---- Dump the restart file when required ------------------------------
*=======================================================================

      IF(nsave .NE. 0) THEN
         IF(mstep .NE. 0 .AND. MOD(mstep,nsave).EQ.0)THEN
            grflag= igr .OR. gofr
            write(kprint,70200)
            CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co,nprot
     &           ,protl,mass,tmass,vcax,vcay,vcaz)

            IF(save_all_restarts) THEN
              nrestart=nrestart+1
              write(chari,'(i4.4)') nrestart 
#ifdef _MPI_
              ibatt=int((para_index-1)/ntraj_eff) + 1
              write(bchari,'(i2.2)') ibatt
              if(nbatteries.gt.1) THEN 
                name_rest= TRIM(restart_out)//bchari//'_'//chari//'.rst'
              ELSE
                name_rest= TRIM(restart_out)//chari//'.rst'
              endif
#else
              name_rest= TRIM(restart_out)//chari//'.rst'
#endif
              CALL dumprs(kdump_out,name_rest,nstep,temp
     &             ,ntap,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,xpg,ypg
     &             ,zpg,xpcma,ypcma,zpcma,vpx,vpy,vpz,vcax,vcay,vcaz,eta
     &             ,etap,vh1,3,thermos,sumarray,ssmarray,navg,navg,co,oc
     &             ,cpress,vco,zz1,lzz,lzz1,grflag,krdf,igrn
     &             ,save_all_restarts)
            ELSE
              CALL dumprs(kdump_out,restart_out,nstep
     &             ,temp,ntap,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,xpg
     &             ,ypg,zpg,xpcma,ypcma,zpcma,vpx,vpy,vpz,vcax,vcay,vcaz
     &             ,eta,etap,vh1,3,thermos,sumarray,ssmarray,navg,navg
     &             ,co,oc,cpress,vco,zz1,lzz,lzz1,grflag,krdf,igrn
     &             ,save_all_restarts)
            END IF
         END IF
      END IF

*=======================================================================
*--- Stop the run smoothly if a file called STOP is found --------------
*=======================================================================

      INQUIRE(FILE='./STOP',EXIST=exist)
      IF(exist) THEN
         grflag= igr .OR. gofr   
         write(kprint,70200)
         CALL matinv(3,3,co,oc,volume)
         volume=volume*boxl**3
         CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co,nprot
     &        ,protl,mass,tmass,vcax,vcay,vcaz)
         CALL dumprs(kdump_out,restart_out,nstep,temp
     &        ,ntap,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,xpg,ypg,zpg
     &        ,xpcma,ypcma,zpcma,vpx,vpy,vpz,vcax,vcay,vcaz,eta,etap,vh1
     &        ,3,thermos,sumarray,ssmarray,navg,navg,co,oc,cpress,vco
     &        ,zz1,lzz,lzz1,grflag,krdf,igrn,save_all_restarts)
         IF(wrtgyr) THEN
            CLOSE(kgyr)
         END IF
#ifdef _MPI_
         call MPI_FINALIZE(ierr)
#endif
         write(kprint,70700)
         STOP
      END IF
 
*=======================================================================
*---- Print out all energies and forces if (ltest_times)            ----
*=======================================================================

      if(nplot_steer.gt.0) THEN 
        naus=mod(mstep,nplot_steer)  
      ELSE 
        naus= 1
      END IF
#ifdef _MPI_ 
      IF( ltest_times .or. nplot_steer.gt.0 .or. rem_run .or.
     &     sim_tempering) THEN
#else
        IF( ltest_times .or. nplot_steer.gt.0 .or. sim_tempering) THEN
#endif         
        CALL kinetic(ss_index,co,nato_slv,nmol,cnstpp_slv,ntap
     &        ,cnstpp,vpx,vpy,vpz,vcax,vcay,vcaz,nprot,protl,mass
     &       ,tmass,wmtp,ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm
     &       ,stresstk,cpress,isostress,vco,masspp,ucepr,temppra
     &       ,massinfty)
         IF(thermos) THEN
            CALL comp_thermos_energy(neta,ndf_thermos,t,qmass,eta
     &           ,etap,uceh,hpot,temph)
         END IF

         call steer_variables(nbonds_added,nbends_added
     &        ,nitors_added,atom_b1,atom_b2,atom_be1,atom_be2,atom_be3
     &        ,atom_it1,atom_it2,atom_it3,atom_it4,xp0,yp0,zp0
     &        ,bond_added,bend_added,tors_added,potit(1,2),litor
     &        ,nitors_added)
         CALL mts_test(efact,ntmtss,fpx_h,fpx_l,fpy_m,fpx_n1,ntmtsp,time
     &        ,urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m,ucns_h,ucns_l
     &        ,ucns_m,uslvbon,uslvben,uslvtor,uslvitor ,conf_bnd_slv_n1
     &        ,coul_bnd_slv,coul_bnd_slv_n1,self_slv ,ucek,urcp_h,urcp_l
     &        ,urcp_m,ucop_h,ucop_l,ucop_m,ucnp_h ,ucnp_l,ucnp_m,ubond
     &        ,ubend,uitors,uptors,fsbond,fsbend ,fsin14,cnb14,unb14
     &        ,cngrp,ungrp,conf_bnd_slt_n1 ,coul_bnd_slt,coul_bnd_slt_n1
     &        ,self_slt,pucek,urcsp_h ,urcsp_l,urcsp_m,ucosp_h,ucosp_l
     &        ,ucosp_m,ucnsp_h,ucnsp_l ,ucnsp_m,eer_m,eer_l,eer_h,nmol
     &        ,nato_slv,ntap,nstep,ktest ,lfirst,maxstp,cpress,volume
     &        ,pext,ucepr,thermos,uceh,hpot,eqdist,eqdist1,force_cost
     &        ,strbonds,nbonds_added,eqang,eqang1,force_ang,strbends
     &        ,nbends_added,eqdied,eqdied1,force_died,strtors
     &        ,nitors_added,str_vel,ben_vel,tor_vel,str_time0,str_time1
     &        ,str_tim0r,str_tim1r,rtime,kplot_steer,bond_added
     &        ,bend_added,tors_added,nstep_steer,lstretch,lbend,litor
     &        ,nplot_steer,time,ltest_times,vi,chrge,steer_com,workcom
     &        ,bussi,ebussi,itprint,ifprint)
      END IF
