      SUBROUTINE mts_furier(xpb,ypb,zpb,xp0,yp0,zp0,xpcm,ypcm,zpcm,urcsp
     &     ,urcs,urcp,virsp,virs,virp,fpx,fpy,fpz,fsin14,gsin14,fsbend
     &     ,gsbend,fsbond,gsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt
     &     ,coul_bnd_slv,rsh,rshk,eer,virial,fudgec)

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

*-------------------- FUNCTION DEFINITION ------------------------------

      EXTERNAL near0
      LOGICAL near0

*-------------------- LOCAL STUFF --------------------------------------

      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      REAL*8   virial(3,3),eer
      INTEGER  i,j,nbond_aux,iat

*---- VARIABLES IN INCLUDE --------------------------------------------*

      REAL*8  fx(m1),fy(m1),fz(m1),fxx(m1,2),fyy(m1,2),fzz(m1,2)
      COMMON /rag1/ fx,fy,fz,fxx,fyy,fzz

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
      IF(rsh .NE. rshk) RETURN
      if(omp_timing) call timer(treal,tcpu1,tela1)

      IF(.NOT. pme) THEN
         CALL furipp(ss_index,oc,xp0,yp0,zp0,chrge,ntap,atomp,grppt
     &        ,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm,zpcm,fpx
     &        ,fpy,fpz,co,virial)
      ELSE

*=======================================================================
*--- Call fft_pme                                                   ----
*=======================================================================
      
         call fft_pme(ntap,xp0,yp0,zp0,xpcm,ypcm,zpcm,chrge,co,oc
     &        ,volume,alphal,pme_order,nfft1,nfft2,nfft3,eer,fpx,fpy
     &        ,fpz,virial,atomp,grppt,pressure,rkcut)
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

*=======================================================================
*---- Compute stress tensor if coupling is by group --------------------
*=======================================================================

c=====      subtract "intramolecular term" in the ZERO cell: 
c--         bonds

      if(lstretch.GT.0) CALL ferrf(ss_index,alphal,chrge,1.0D0,xpb,ypb
     &     ,zpb,1,lstrtch,lstretch,fsbond_slt,fsbond_slv,fpx,fpy,fpz
     &     ,erf_corr,erf_arr_corr,delew,rlew,nbonds_added,ntap)
      fsbond=fsbond_slt+fscnstr_slt
      coul_bnd_slv=fsbond_slv+fscnstr_slv
      coul_bnd_slt=fsbond_slt+fscnstr_slt
      
c--         bendings

      nbond_aux=0
      CALL ferrf(ss_index,alphal,chrge,1.d0,xpb,ypb,zpb,1,int13,int13p
     &     ,fsbend_slt,fsbend_slv,fpx,fpy,fpz,erf_corr,erf_arr_corr
     &     ,delew,rlew,nbond_aux,ntap)

      coul_bnd_slt=coul_bnd_slt+fsbend_slt
      coul_bnd_slv=coul_bnd_slv+fsbend_slv
      fsbend=fsbend_slt

c---         torsions (fudged)
      
      CALL ferrf(ss_index,alphal,chrge,fudgec,xpb,ypb,zpb,1,int14,int14p
     &     ,fsin14_slt,fsin14_slv,fpx,fpy,fpz,erf_corr,erf_arr_corr
     &     ,delew,rlew,nbond_aux,ntap)
      fsin14=fsin14_slt
      coul_bnd_slv=coul_bnd_slv+fsin14_slv
      coul_bnd_slt=coul_bnd_slt+fsin14_slt

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_f=time_f+tela2-tela1
        tcpu_f=tcpu_f+tcpu2-tcpu1
      END IF
      RETURN
      END
