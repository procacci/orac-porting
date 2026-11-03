      subroutine mts_intra_n1(xp0,yp0,zp0,xcm,ycm,zcm,fpx,fpy,fpz
     &     ,fudge_qq,fudge_lj,fold_dir,puhyd,conf_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slt,coul_bnd_slv,unb14,cnb14,ungrp,cngrp,uptors
     &     ,uslvtor,mapdn,nmapdn,uumb,gr)

************************************************************************
*     FINTRAP call slow "bonded" force routines for the solute.   
************************************************************************

      use parst
      use cpropar
      use omp_integr
      use rem
#ifdef _OMP_
#ifdef _BGQ_  
      use omp_lib       ! OMP layer for BGQ only
#endif
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  xp0(*),yp0(*),zp0(*),xcm(*),ycm(*),zcm(*),unb14,cnb14
      REAL*8  fpx(ntap),fpy(ntap),fpz(ntap),puhyd
     &     ,conf_bnd_slt,conf_bnd_slv,coul_bnd_slt,coul_bnd_slv,uslvtor
      INTEGER fold_dir,nmapdn(*),mapdn(2,*) 
      REAL*8 uptors,fudge_qq,fudge_lj,ungrp,cngrp,uumb,gr

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  unb14_slt,cnb14_slt,ungrp_slt,cngrp_slt,uptors_slt
     &       ,unb14_slv,cnb14_slv,ungrp_slv,cngrp_slv,uptors_slv
      REAL*8  xc,yc,zc,st(3,3),seg_ene1,seg_ene2
      INTEGER i,j,i1,count,m,iat
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif
      integer itask

      if(omp_timing) THEN 
        call timer(treal,tcpu1,tela1)
      end if

#ifdef _OMP_
      if(omp_dynamic) call omp_set_dynamic(.TRUE.)
      unb14_slt=0.d0
      unb14_slv=0.d0
      cnb14_slt=0.d0
      cnb14_slv=0.d0
      ungrp_slt=0.d0
      ungrp_slv=0.d0
      cngrp_slt=0.d0
      cngrp_slv=0.d0
      uptors_slt=0.d0
      uptors_slv=0.d0
      dwrk=0.d0
!$OMP    PARALLEL  PRIVATE(ITASK) NUM_THREADS(nthr2)
!$OMP&   REDUCTION(+:unb14_slt,unb14_slv,cnb14_slt,cnb14_slv)
!$OMP&   REDUCTION(+:ungrp_slt,ungrp_slv,cngrp_slt,cngrp_slv)
!$OMP&   REDUCTION(+:uptors_slt,uptors_slv,fpx,fpy,fpz)
!$OMP&   REDUCTION(+:seg_ene1,seg_ene2,dwrk)
      itask=1+OMP_GET_THREAD_NUM()
#else
      itask=1
#endif

*=======================================================================
*----- Third neighbour interactions : Sum on the direct lattice --------
*=======================================================================
      if(rem_segment) then
        call fnb14_rem(ss_index,xp0,yp0,zp0,chrge,ntap,ecc1412,ecc146
     &       ,cutoff,clewld,alphal,int14,int14p,type14,fudge,lj_fudge
     &       ,fpx,fpy,fpz,unb14_slt,unb14_slv,cnb14_slt,cnb14_slv
     &       ,remint14,itask,seg_ene1)
      else
        if(alchemy) THEN 
          call fnb14_alchemy(ss_index,xp0,yp0,zp0,chrge,ntap,ecc1412
     &         ,ecc146,cutoff,clewld,alphal,int14,int14p,type14,fudge
     &         ,lj_fudge,fpx,fpy,fpz,unb14_slt,unb14_slv,cnb14_slt
     &         ,cnb14_slv,lambda,lambdaq,lambda0,lambdaq0,dwrk,itask)
        ELSE
          call fnb14(ss_index,xp0,yp0,zp0,chrge,ntap,ecc1412
     &         ,ecc146,cutoff,clewld,alphal,int14,int14p,type14,fudge
     &       ,lj_fudge,fpx,fpy,fpz,unb14_slt,unb14_slv,cnb14_slt
     &         ,cnb14_slv,itask)
        END IF
      endif

         
*=======================================================================
*----- Interactions between atoms contained in the same group ----------
*=======================================================================

      if(alchemy) THEN 
        CALL fnbgrp_alchemy(ss_index,xp0,yp0,zp0,chrge,nbtype,ecc12,ecc6
     &       ,clewld,alphal,ingrp,ingrpp,fpx,fpy,fpz,ungrp_slt,ungrp_slv
     &       ,cngrp_slt,cngrp_slv,lambdaq,itask)
      ELSE
        CALL fnbgrp(ss_index,xp0,yp0,zp0,chrge,nbtype,ecc12,ecc6
     &       ,clewld,alphal,ingrp,ingrpp,fpx,fpy,fpz,ungrp_slt,ungrp_slv
     &       ,cngrp_slt,cngrp_slv,ntap,itask)
      END IF


*=======================================================================
*----- Bonded interactions: proper torsions ----------------------------
*=======================================================================      
      IF(rem_segment) THEN
        CALL fptors_rem(ss_index,ltor,ltors,xp0,yp0,zp0,potto(1,2)
     &       ,potto(1,1),uptors_slt,uptors_slv,fpx,fpy,fpz,remltor
     &       ,itask,seg_ene2)
      ELSE
        CALL fptors(ss_index,ltor,ltors,xp0,yp0,zp0,potto(1,2),potto(1
     &       ,1),uptors_slt,uptors_slv,fpx,fpy,fpz,m1,itask)
      ENDIF
#ifdef _OMP_
!$OMP   END PARALLEL
      if(omp_dynamic) THEN
        call omp_set_dynamic(.FALSE.)
        call OMP_SET_NUM_THREADS(nthr) !restore thread numebr
      END IF
#endif
      if(rem_segment) rem_n1_energy=seg_ene1 + seg_ene2
      unb14=unb14_slt
      cnb14=cnb14_slt
      
      conf_bnd_slt=unb14_slt+ungrp_slt
      coul_bnd_slt=cnb14_slt+cngrp_slt
      conf_bnd_slv=unb14_slv+ungrp_slv
      coul_bnd_slv=cnb14_slv+cngrp_slv

      ungrp=ungrp_slt
      cngrp=cngrp_slt

      uptors =uptors_slt
      uslvtor=uptors_slv

*=======================================================================
*---- Add forces from hystory-dependent potential if present ----
*=======================================================================
      
      call meta_forces(xp0,yp0,zp0,fpx,fpy,fpz,co,oc,ntap)
      
c--   set to zero forces on fixed atoms
      do i=1,nfixed
        iat=ifixed(i)
        fpx(iat)=0.d0
        fpy(iat)=0.d0
        fpz(iat)=0.d0
      end do
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2) 
        time_n1=time_n1+tela2-tela1
        tcpu_n1=tcpu_n1+tcpu2-tcpu1
      end if
      RETURN
      END
