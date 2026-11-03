      SUBROUTINE get_total_energy(i_forces,mapnl,mapdn,nmapdn,nnlpp0
     &     ,nnlpp,abmd_dir,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,prt
     &     ,utotal,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,self_slv,fsin14,unb14,cnb14,uslvbon,uslvben
     &     ,uslvtor,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop
     &     ,urcp,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt
     &     ,ucnsp,ucosp,urcsp,eer,enisolv,nstep)

************************************************************************
*   Time-stamp: <98/02/21 11:05:29 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 14 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar
      use sasamod
#ifdef _OMP_
      use omp_integr
#endif

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER abmd_dir,mapdn(2,*),nmapdn(*)
#ifdef _OMP_
      INTEGER  nnlpp0(mpp8,nthr),nnlpp(mpp8,nthr),mapnl(m8t,nthr)
     &     ,ncount1(nthr)
#else
      INTEGER nnlpp0(*),nnlpp(*),mapnl(*)
#endif
      INTEGER nstep
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
     &     ,prt(3,3),utotal
      LOGICAL i_forces

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8 xpa(m1),ypa(m1),zpa(m1),xpcma(npm),ypcma(npm),zpcma(npm)
     &     ,xpcm(npm),ypcm(npm),zpcm(npm),xpga(m11),ypga(m11),zpga(m11)
     &     ,xpg(m11),ypg(m11),zpg(m11),fpx1(m1),fpy1(m1),fpz1(m1)
     &     ,prt1(3,3)

      COMMON /rag2/ fpx1,fpy1,fpz1,prt1,xpa,ypa,zpa,xpcm,ypcm,zpcm,xpcma
     &     ,ypcma,zpcma,xpga,ypga,zpga,xpg,ypg,zpg

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors,uconf,ucoul
     &     ,ureal,urecp,urcs,urcp,urcsp,urcsp_m,eer,ucns,ucos,ucnsp
     &     ,ucosp,ucnp,ucop,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,virs,virsp,virp,gsin14,fsbend,fsbond,gsbend
     &     ,gsbond,unb14,cnb14,fsin14,ungrp,cngrp,uumb,gr,ubond,uslvbon
     &     ,uslvben,unbond,cnbond,purecp,upconf,upcoul,fudgec,press
     &     ,stressd(3,3),stressr(3,3),st(3,3),pv,ustot,uptot,upstot
     &     ,enisolv,aux
      INTEGER iret,i
      CHARACTER*80 errmsg

      CHARACTER*1 rshell,rshk
      DATA   rshell/'h'/rshk/'h'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot
     &     ,protl)

*=======================================================================
*---- Calculate group position  ----------------------------------------
*=======================================================================

      CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
      
*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================

      CALL fndgrp(nprot,protl,atomp)

*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &     ,zpcma)

      urcs=0.0D0
      urcp=0.0D0
      urcsp=0.0D0
      eer=0.0D0
      self_slt=0.0D0
      self_slv=0.0D0
      fscnstr_slt=0.0D0
      fscnstr_slv=0.0D0
      fsin14=0.0D0
      unb14=0.0D0
      cnb14=0.0D0
      coul_bnd_slt=0.0D0
      coul_bnd_slv=0.0D0
      CALL zeroa(fpx1,fpy1,fpz1,ntap,1)
      IF(clewld) THEN
         CALL cself(ss_index,ntap,alphal,rkcut,chrge,self_slt,self_slv)
         CALL ferrf(ss_index,alphal,chrge,1.0D0,xp0,yp0,zp0,0,lcnstr
     &        ,lconstr,fscnstr_slt,fscnstr_slv,fpx1,fpy1,fpz1
     &        ,erf_corr,erf_arr_corr,delew,rlew,nbonds_added,ntap)
      END IF

      CALL mts_forces(nstep,rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns,ucos,virs,virsp,ucnp
     &     ,ucop,ucnsp,ucosp,fpx1,fpy1,fpz1,stressd,nnlpp0,nnlpp)

      IF(clewld) THEN
         CALL mts_furier(xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp
     &        ,urcs,urcp,virsp,virs,virp,fpx1,fpy1,fpz1,fsin14
     &        ,gsin14,fsbend,gsbend,fsbond,gsbond,fscnstr_slt
     &        ,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell,rshk,eer
     &        ,stressr,fudgec)
      END IF



      CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx1,fpy1
     &     ,fpz1,fudge,lj_fudge,abmd_dir,puhyd,conf_bnd_slt_n1
     &     ,conf_bnd_slv_n1,coul_bnd_slt_n1,coul_bnd_slv_n1,unb14,cnb14
     &     ,ungrp,cngrp,uptors,uslvtor,mapdn,nmapdn,uumb,gr)
      CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx1,fpy1
     &     ,fpz1,ubond,uslvbon,ubend,uslvben,uitors,uslvitor
     &     ,ntap)

      IF(slv_exist) THEN
         uconf=ucns+ucos+urcs+coul_bnd_slv+conf_bnd_slv_n1
     &        +coul_bnd_slv_n1+self_slv + uumb
         ucoul=ucos+urcs+coul_bnd_slv+coul_bnd_slv_n1+self_slv
         ureal=ucos+coul_bnd_slv_n1
         urecp=urcs+self_slv+coul_bnd_slv
      END IF
      
*---        Energy terms for the protein         -----------------------
      
      IF(slt_exist) THEN 
         pubnd =uptors+uitors+ubond+ubend
         unbond=ucnp
         cnbond=ucop
         purecp=urcp
         puconf= unbond + conf_bnd_slt_n1 + uumb
         pucoul= cnbond + purecp + coul_bnd_slt + coul_bnd_slt_n1
     &        +self_slt
      END IF   
      
      
*---        Mixed terms                          -----------------------
      
      IF(slv_exist .AND. slt_exist) THEN 
         upconf= ucnsp+ucosp+urcsp
         upcoul= ucosp+urcsp
      END IF
      
      IF(pressure) THEN
         CALL comp_stress_conf(stressd,stressr,prt,oc
     &        ,volume,unitp,press)
         CALL zero0(st,9)
         CALL comp_forcep(prt,st,oc,volume,pext)
      END IF
      
*----       If pme=T no solute-solvent separation ----------------------
      
      IF(pme) THEN 
         IF(slv_exist.AND.(.NOT.slt_exist)) THEN
            uconf=uconf + eer
            ucoul=ucoul + eer
            urecp= eer 
         END IF
         IF(slt_exist.AND.(.NOT.slv_exist)) THEN
            purecp=eer 
            pucoul=pucoul+eer 
         END IF
         IF(slt_exist.AND.slv_exist) THEN 
            upconf=upconf + eer 
            upcoul=upcoul + eer 
         END IF
      END IF
 
      pv=0.0D0
      IF(pressure .OR. cpress .OR. lberendsen) THEN
         pv=(pext*volume*efact)/1000.0D0
      END IF
      ustot=0.d0
      ustot=(ucns+ucos+urcs+uslvbon+uslvben+uslvtor+uslvitor
     &     +conf_bnd_slv_n1+coul_bnd_slv+coul_bnd_slv_n1+self_slv)*efact
     &     /1000.0D0
      uptot=0.d0
      uptot=(urcp+ucop+ucnp+ubond+ubend+uptors+uitors
     &     +conf_bnd_slt_n1+coul_bnd_slt+coul_bnd_slt_n1+self_slt)*efact
     &     /1000.d0

      upstot=0.d0
      upstot = (urcsp+ucosp
     &     +ucnsp)*efact/1000.d0

      utotal= ustot+uptot+upstot+eer*efact/1000.d0
      IF(cpress) THEN
         pv=(pext*volume*efact)/1000.0D0
         utotal=utotal+pv
      END IF
      if(agbnp) THEN 
        aux = 4.814d0*1000.d0/efact 
        call agbnpf(xp0,yp0,zp0,fpx1,fpy1,fpz1,aux,enisolv,iret) ! call agbnp interface
        utotal=utotal+enisolv*efact/1000.d0
      else if(lsasa) THEN 
        call force_sasa(ntap,sasa_nbtype,sasa_betb,xpa,ypa,zpa,co,r2cut
     &       ,rsolv,rvdw,ssasa,psasa,enisolv,fpx1,fpy1,fpz1,sasa)
        utotal=utotal+enisolv*efact/1000.d0
      END IF
      IF(i_forces) THEN
         CALL dcopy(ntap,fpx1,1,fpx,1)
         CALL dcopy(ntap,fpy1,1,fpy,1)
         CALL dcopy(ntap,fpz1,1,fpz,1)
      END IF 

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
