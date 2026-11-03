!===============================================================================================================================
      SUBROUTINE  mts_forces(nstep,rshell,xp0,yp0,zp0,xpg,ypg,zpg,xpcm
     &     ,ypcm,zpcm,mapnl,mapdn,nmapdn,ucns,ucos,virs,virsp,ucnp,ucop
     &     ,ucnsp,ucosp,fpx,fpy,fpz,stress,nnlppo,nnlppi)
!===============================================================================================================================
!                 written by P Procacci Florence Univ 2014                                                          
!     Transit force routine: calls all non bonded force routines (REM Segment, Alchemy, standard)
!===============================================================================================================================

      use unit
      use parst
      use cpropar
      use rem
      use omp_integr
      
  !======================= DECLARATIONS ==================================

      IMPLICIT none

  !----------------------- ARGUMENTS -------------------------------------

      CHARACTER*1  rshell
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xpg(*),ypg(*)
     &     ,zpg(*),xpcm(*),ypcm(*),zpcm(*),stress(3,3)
      REAL*8  ucnp,ucop,ucns,ucos,ucnsp,ucosp,virs,virsp

#ifdef _OMP_
      INTEGER mapnl(m8t,nthr),mapdn(2,*),nmapdn(*),nstep,nnlppo(mpp8
     &     ,nthr),nnlppi(mpp8,nthr),i
#else      
      INTEGER mapnl(*),mapdn(2,*),nmapdn(*),nstep,nnlppo(*),nnlppi(*)
#endif

  !-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      EXTERNAL  near0
      LOGICAL near0

  !-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 ro,rto,ri,rti,rni,seg_ene
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8 treal
#ifdef _OMP_
      INTEGER npp(nthr)
      real*8 nvtot_av,dnvtot_av
#else
      INTEGER npp
#endif
      INTEGER iz1,m,n,naux
      LOGICAL lskip_ewald
      
  !==================== EXECUTABLE STATEMENTS ============================

  ! select cutoffs according to current value of SHELL
      
      if(rshell.eq.'u') then
        if(omp_timing) call timer(treal,tcpu1,tela1)
        iz1=0
        ro=rcuth
        ri=rcuth
        rto=0.d0
        rti=rtolh
        rni=rneih
      endif
      if(rshell.eq.'h') then
        if(omp_timing) call timer(treal,tcpu1,tela1)
        iz1=1
        if(rneih.lt.1.d-10) iz1=0
        ro=rcuth
        ri=rcutl
        rto=rtolh
        rti=rtoll
        rni=rneil
      endif
      if(rshell.eq.'l') then
        if(omp_timing) call timer(treal,tcpu1,tela1)
        iz1=1
        ro=rcutl
        rto=rtoll
        ri=rcutm
        rti=rtolm
        rni=rneim
      end if
      if(rshell.eq.'m') then
        if(omp_timing) call timer(treal,tcpu1,tela1)
        iz1=1
        ro=rcutm
        rto=rtolm
        ri=0.d0
        rti=0.d0
        rni=0.d0
      end if
      
!     set energies to zero 
      ucnp=0.d0
      ucop=0.d0
      ucns=0.d0
      ucos=0.d0
      virs=0.d0
      ucnsp=0.d0
      ucosp=0.d0
      DO m=1,3
        DO n=1,3
          stress(n,m)=0.0D0
        END DO
      END DO
      
! CALL to Force routines
          
      if(clewld) THEN 
        lskip_ewald=.false.   
        call decide_force(ri,alphal,lskip_ewald)
      ELSE
        lskip_ewald=.true.
      END IF
  
      IF(rem_segment) THEN      !====== REM (segment) =========
! only lj interactions in long-range/medimum shell
        call force_rem(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,chrge,nbtype
     &       ,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &       ,ecc12,ecc6,alphal,mapnl,ngrp,grppt,ucnp,ucns,ucnsp,ucop
     &       ,ucos,ucosp,fpx,fpy,fpz,stress,nnlppo,nnlppi,rni,ri,ro,rti
     &       ,rto,pmass,iz1,npp,mass,massinfty,seg_ene,rematom
     &       ,lskip_ewald,clewld) 
        if(rshell.eq.'h') rem_hnb_energy = seg_ene
        if(rshell.eq.'l') rem_lnb_energy = seg_ene
        if(rshell.eq.'m') rem_mnb_energy = seg_ene
      ELSE IF (alchemy) THEN    ! ======== ALCHEMY ===========    
        IF(shifted) THEN        ! --- shifted potential
          call force_ALCHEMY_shift(ss_index,xp0,yp0,zp0,xpg,ypg,zpg
     &         ,chrge,nbtype,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm
     &         ,groupp,atomp,co,ecc12,ecc6,alphal,mapnl,ngrp,grppt,ucnp
     &         ,ucns,ucnsp,ucop,ucos,ucosp,fpx,fpy,fpz,stress,nnlppo
     &         ,nnlppi,rni,ri,ro,rti,rto,pmass,iz1,npp,mass,massinfty
     &         ,lambda,lambda0,lambdaq,lambdaq0,rmin,epsm,dwrk
     &         ,lskip_ewald,clewld,qshift,ljshift)
        ELSE              ! soft-core
          call force_ALCHEMY(ss_index,xp0,yp0,zp0,xpg,ypg,zpg
     &         ,chrge,nbtype,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm
     &         ,groupp,atomp,co,ecc12,ecc6,alphal,mapnl,ngrp,grppt,ucnp
     &         ,ucns,ucnsp,ucop,ucos,ucosp,fpx,fpy,fpz,stress,nnlppo
     &         ,nnlppi,rni,ri,ro,rti,rto,pmass,iz1,npp,mass,massinfty
     &         ,lambda,lambda0,lambdaq,lambdaq0,rmin,epsm,dwrk
     &         ,lskip_ewald,clewld,qsoft,ljsoft)
        END IF
      ELSE                      ! ======= STANDARD ===========  
        call force(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,chrge,nbtype
     &       ,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &       ,ecc12,ecc6,alphal,mapnl,ngrp,grppt,ucnp,ucns,ucnsp,ucop
     &       ,ucos,ucosp,fpx,fpy,fpz,stress,nnlppo,nnlppi,rni,ri,ro,rti
     &       ,rto,pmass,iz1,npp,mass,massinfty,lskip_ewald,clewld)
      ENDIF

      if(nupdte*lrespa*mrespa.gt.nprint) THEN 
#ifdef _OMP_
        if(iz1.eq.0) THEN 
          WRITE(kprint,*)
          call  ave_neigh(nthr,nvtot_av,dnvtot_av,npp)
          IF(nthr.GT.1) THEN
            WRITE(kprint,10001) nvtot_av,dnvtot_av
          ELSE
            WRITE(kprint,10000) nint(nvtot_av)
          END IF
        END IF
#else      
        IF(iz1 .EQ. 0) THEN
          WRITE(kprint,10000) npp
        END IF
#endif
      else 
#ifdef _OMP_
        naux=nupdte*int((nprint/lrespa/mrespa)/nupdte)
        if(mod(nstep,naux).eq.0.and.iz1.eq.0) THEN 
          WRITE(kprint,*)
          call  ave_neigh(nthr,nvtot_av,dnvtot_av,npp)
          IF(nthr.GT.1) THEN
            WRITE(kprint,10001) nvtot_av,dnvtot_av
          ELSE
            WRITE(kprint,10000) nint(nvtot_av)
          END IF
        END IF
#else
        naux=nupdte*int((nprint/lrespa/mrespa)/nupdte)
        if(mod(nstep,naux).eq.0.and.iz1.eq.0) WRITE(kprint,10000)
     &       npp
#endif
      end if

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        if(rshell.eq.'u') time_u=time_u+tela2-tela1
        if(rshell.eq.'h') time_h=time_h+tela2-tela1
        if(rshell.eq.'l') time_l=time_l+tela2-tela1
        if(rshell.eq.'m') time_m=time_m+tela2-tela1
        if(rshell.eq.'u') tcpu_u=tcpu_u+tcpu2-tcpu1
        if(rshell.eq.'h') tcpu_h=tcpu_h+tcpu2-tcpu1
        if(rshell.eq.'l') tcpu_l=tcpu_l+tcpu2-tcpu1
        if(rshell.eq.'m') tcpu_m=tcpu_m+tcpu2-tcpu1
      END IF
      
!================= END OF EXECUTABLE STATEMENTS ========================
      
10000 FORMAT(5x,'Neighbor Lists Dimensions     *neighbor(',i7,')* ')
10001 FORMAT(5x,'Average Neighbor Lists Dimension -> ',F15.1,' +/-'
     &     ,F15.1)

      RETURN
      END       
       
!===========================================================
      subroutine decide_force(rinn,alphal,lskip_ewald)
!===============================================================================================================================
!      decide whether Coulomb forces should be computed according to outer cut-off.
!===============================================================================================================================
!======================= DECLARATIONS ==================================
  
      IMPLICIT none
  
!----------------------- ARGUMENTS -------------------------------------
      REAL*8 rinn,alphal
      LOGICAL lskip_ewald
      
!-!-------------------- LOCAL VARIABLES ----------------------------------
      REAL*8  a1,a2,a3,a4,a5,qp,eps1,erftbdns,qt,expcst,erfcst,rsp,c1,c2
     &     ,c3,c4,h
      INTEGER ic
      REAL*8 zero,one,two,twoi,three,threei
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,      
     &     three=3.0D0,twoi=0.5D0,threei=1.0D0/3.0D0)
      
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

      qt=1.0d0/(1.0d0+qp*alphal*rinn)
      expcst=DEXP(-alphal*rinn*alphal*rinn)
      erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
      lskip_ewald = erfcst.lt.10.d-4
      RETURN
      END 
      
      subroutine ave_neigh(nthr,nvtot_av,dnvtot_av,npp)

      implicit none
      real*8  nvtot_av,dnvtot_av,nvtot_av2 
      integer npp(*),i,nthr

      nvtot_av=0.d0
      nvtot_av2=0.d0
      DO I=1,nthr
        nvtot_av=nvtot_av+npp(i)
        nvtot_av2=nvtot_av2+npp(i)**2.
      end do
      nvtot_av=nvtot_av/float(nthr)
      nvtot_av2= nvtot_av2/float(nthr)
      dnvtot_av=dsqrt(nvtot_av2-nvtot_av**2)
      return
      end 
