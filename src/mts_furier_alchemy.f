      SUBROUTINE mts_furier_alchemy(xpb,ypb,zpb,xp0,yp0,zp0,xpcm,ypcm
     &     ,zpcm,urcsp,urcs,urcp,virsp,virs,virp,fpx,fpy,fpz,fsin14
     &     ,gsin14,fsbend,gsbend,fsbond,gsbond,fscnstr_slt,fscnstr_slv
     &     ,coul_bnd_slt,coul_bnd_slv,rsh,rshk,eer,virial,fudgec
     &     ,self_slt)

************************************************************************
*                                                                      *
*     Transit routine : calls MTS_FURIWW, MTS_FURIPP, MTS_FURIPW,      *
*     FERFF, WATSELF, PME_FFT                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*                                                                      *
*                                                                      *
************************************************************************
! This is the alchemic version. let Q and q be the full and alchemic
! solute charges and S be the full "solvent" charges (that may include
! part of the non alchemic solute)
! In the direct and Recp lattice the following interactions are evaluated:
!     Direct Lattice                           Reciprocal lattice
!   QQ*erfc   qS*erfc   SS*erfc           qq*erf    qS*erf    SS*erf   
!     only 1-4f and >14                            ALL
! Then we need to take away qq erf(12 13 and 14(1-f) for the recp.
! lattice. This is done here, in this very subroutine. Then we need to
! *add* (QQ-qq)*erf in the direct lattice in the shell m for all
! intrasolute such that ij >14. The fudged term is take into care in 
! fnb14_alchemy. The qq self term is automatically included in the recp
! and must be also removed. This is done here. 
! The term that do work are:
! Direct lattice:  qS [only in mts_alchemy_forpp]
! Reciprocal Lattice: pme, -ferrf(12) ,-ferrf(13), -ferrf(14)(1-f) -self
! [HERE]  and -ferrf(>14) in mtsmd.f and -ferrf*f in fnb14_alchemy.
! NB for the intra-alchemical solute ALL the 14 qq interactions must be
! subtracted as they must be QQ*f/r (i.e. [erfc +erf]QQ/r. Part is
! done here (1-f) and part in the fnb14_alchemy (f term only)    

*======================= DECLARATIONS ==================================

      use parst
      use cpropar
      use omp_integr
      use spme

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8     xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xpb(*),ypb(*)
     &     ,zpb(*)
      REAL*8     urcsp,urcs,urcp,virsp,virs,virp,fsin14,gsin14,fsbend
     &     ,gsbend,fsbond,gsbond,fudgec,fscnstr_slt,fscnstr_slv
     &     ,coul_bnd_slt,coul_bnd_slv
      character*1 rsh,rshk
      REAL*8  xpcm(*),ypcm(*),zpcm(*)
      REAL*8   virial(3,3),eer,self_slt

*-------------------- FUNCTION DEFINITION ------------------------------

      EXTERNAL near0
      LOGICAL near0

c---  stuff for b-spline interpolation and FFT 

c---  local stuff 
      
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      INTEGER  i,j,nbond_aux,iat
      REAL*8 chrge0(m1),chrgec(m1),dwrk_bo,dwrk_cs,dwrk_be,dwrk_to
     &     ,dwrk_al,dwrk_rl,dwrk_sf,fx0,fy0,fz0

*---- VARIABLES IN INCLUDE --------------------------------------------*

      REAL*8  fsbond_slt,fsbend_slt,fsin14_slt,fsbond_slv,fsbend_slv
     &     ,fsin14_slv

*==== EXECUTABLE STATEMENTS: ==========================================*


c---------------------------------------------------------------------
c     PROTEIN alone or PROTEIN + SOLVENT
c---------------------------------------------------------------------

      DO i=1,3
         DO j=1,3
            virial(i,j)=0.0D0
         END DO
      END DO
      eer=0.0D0
      dwrk=0.d0
      IF(rsh .NE. rshk) RETURN
      if(omp_timing) call timer(treal,tcpu1,tela1)
      coul_bnd_slt=0.d0

c---  set current and previous charges 
      do i=1,ntap 
        if(dabs(lambdaq(i)).eq.1.d0) THEN 
          chrgec(i)=0.d0 
        ELSE if(dabs(lambdaq(i)).eq.0.d0) THEN
          chrgec(i)=chrge(i)
        ELSE
          chrgec(i)=( 1.-dabs(lambdaq(i))) *chrge(i)
        end if
        if(dabs(lambdaq0(i)).eq.1.d0) THEN 
          chrge0(i)=0.d0 
        ELSE if(dabs(lambdaq0(i)).eq.0.d0) THEN
          chrge0(i)=chrge(i)
        ELSE
          chrge0(i)=(1.-dabs(lambdaq0(i)))*chrge(i)
        end if

      end do
      IF(.NOT. pme) THEN
        CALL furipp_alchemy(ss_index,oc,xp0,yp0,zp0,chrgec,chrge0,ntap
     &       ,atomp,grppt,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm
     &       ,zpcm,fpx,fpy,fpz,co,virial,dwrk_rl)
      ELSE

*=======================================================================
*--- Call fft_pme                                                   ----
*=======================================================================
      
        call fft_pme_alchemy(ntap,xp0,yp0,zp0,xpcm,ypcm,zpcm,chrgec
     &       ,chrge0,co,oc,volume,alphal,pme_order,nfft1,nfft2,nfft3,eer
     &       ,fpx,fpy,fpz,virial,atomp,grppt,pressure,rkcut,dwrk_rl)
      END IF 

c--   set to zero forces on fixed atoms
      do i=1,nfixed
        iat=ifixed(i)
        fpx(iat)=0.d0
        fpy(iat)=0.d0
        fpz(iat)=0.d0
      end do
 
      IF(remove_momentum) THEN
         CALL remove_mv(fpx,fpy,fpz,mass,ntap,massinfty)
      END IF

c=====      subtract "intramolecular term" in the ZERO cell: 
c--         bonds + constraints

c---------------------------------------------------------------------------
      CALL ferrf_alchemy(ss_index,alphal,chrgec,chrge0,1.0D0,xpb,ypb,zpb
     &     ,1,lstrtch,lstretch,fsbond_slt,fsbond_slv,fpx,fpy,fpz
     &     ,erf_corr,erf_arr_corr,delew,rlew,nbonds_added,dwrk_bo,ntap)
      CALL ferrf_alchemy(ss_index,alphal,chrgec,chrge0,1.0D0,xpb,ypb
     &     ,zpb,0,lcnstr,lconstr,fscnstr_slt,fscnstr_slv,fx0,fy0
     &     ,fz0,erf_corr,erf_arr_corr,delew,rlew,0,dwrk_cs,ntap)
      fsbond=fsbond_slt+fscnstr_slt

      coul_bnd_slv=fsbond_slv+fscnstr_slv
      coul_bnd_slt=fsbond_slt+fscnstr_slt
c---------------------------------------------------------------------------
      
c--         bendings
c---------------------------------------------------------------------------
      nbond_aux=0
      CALL ferrf_alchemy(ss_index,alphal,chrgec,chrge0,1.D0,xpb,ypb,zpb
     &     ,1,int13,int13p,fsbend_slt,fsbend_slv,fpx,fpy,fpz,erf_corr
     &     ,erf_arr_corr,delew,rlew,nbond_aux,dwrk_be,ntap)
      coul_bnd_slt=coul_bnd_slt+fsbend_slt
      coul_bnd_slv=coul_bnd_slv+fsbend_slv
      fsbend=fsbend_slt
c---------------------------------------------------------------------------

c---------------------------------------------------------------------------
c---         torsions (fudged) 
      nbond_aux=0
      CALL ferrf_alchemy(ss_index,alphal,chrgec,chrge0,fudgec,xpb,ypb
     &     ,zpb,1,int14,int14p,fsin14_slt,fsin14_slv,fpx,fpy,fpz
     &     ,erf_corr,erf_arr_corr,delew,rlew,nbond_aux,dwrk_to,ntap)
      fsin14=fsin14_slt
      coul_bnd_slv=coul_bnd_slv+fsin14_slv
      coul_bnd_slt=coul_bnd_slt+fsin14_slt
c---------------------------------------------------------------------------

c---------------------------------------------------------------------------
c---    qq  self term must be taken away 
        self_slt=0.d0
        CALL cself_alchemy(ss_index,ntap,alphal,rkcut,chrgec,chrge0
     &       ,self_slt,dwrk_sf)

      dwrk=dwrk_bo+dwrk_cs+dwrk_be+dwrk_to+dwrk_rl+dwrk_sf

c#ifdef DEBUGALCHEMY
c      write(96,102) coul_bnd_slt*1.d4,fsin14_slt*1.d4,fsbend_slt*1.d4
c     &     ,fsbond_slt*1.d4,fscnstr_slt*1.d4,self_slt*1.d4
c102   format("COUL_BND_SLT",F12.5," 14",f12.5, " be",f12.5, " bo",f12.5
c     &     , " cs",f12.5," sf",f12.5)
c      write(97,103) dwrk_bo*1.d4,dwrk_cs*1.d4,dwrk_be*1.d4,dwrk_to*1.d4
c     &     ,dwrk_rl*1.d4,dwrk_sf*1.d4
c103   format("WRKS"," bo",G15.5" cs",G15.5" be",G15.5" to",G15.5,
c     &     " rl",G15.5, " sf",G15.5 )
cc#endif
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_f=time_f+tela2-tela1
        tcpu_f=tcpu_f+tcpu2-tcpu1
      END IF
      RETURN
      END
