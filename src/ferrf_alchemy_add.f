      SUBROUTINE ferrf_alchemy_add(ss_index,alphal,charge,lambdaq
     &     ,lambdaq0,x0,y0,z0,list,nlist,fscnstr_slt,fscnstr_slv,fpx,fpy
     &     ,fpz,erf_corr,erf_arr_corr,delew,rlew,dwrk,nato,intramol)

************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald correction deriving from intra alchemical      *
*     interaction for all ij > 14                                      *
*                                                                      *
*   Let Q and q be the full and changing charges of the alchemical     *
*   unit. During the alchemical transformation the Coulomb             *
*   intra alachemical solute(s) must remain unchanged. Only the        *
*   the interaction (alchemical)solute-solvent changes during the      * 
*   charging or  uncharging process if intramol=no (default)                                 * 
*   The qq 12 13 14(1-fudge) interactions are subtracted in the        *
*   mts_furier_alchemy routine. Here, the term (erf)qq/r is subtracted *
*   for >14 contacts  and the erf QQ/r is added to match the bare      *
*   potential when summed to the direct lattice erfc QQ/r computed in  *
*   mts_alchemy_forpp. Due care must be taken to add only intra-solute *
*   terms  erfc QQ/r WITHIN added (or removed) species since added     *
*   species do not free removed species.                               *
*                                                                      *
*---- Last update 27/01/13 --------------------------------------------*
*                                                                      *
*     Written Piero Procacci Chem. dep. 2013                           *
*                                                                      *
*     EXTERNALS  NONE                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
#ifdef _OMP_
      use omp_integr, only:nthr1,nthr,omp_dynamic
#endif
#ifdef _OMP_
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iz,nlist,list(3,*),ss_index(*),nato
      REAL*8  charge(*),lambdaq(*),lambdaq0(*),x0(*),y0(*),z0(*)
     &     ,fpx(nato),fpy(nato),fpz(nato),fscnstr_slt,fscnstr_slv
     &     ,erf_arr_corr(4,*),delew,rlew,dwrk
      REAL*8  alphal
      LOGICAL erf_corr,intramol
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif

*-------------------- LOCAL VARIABLES ----------------------------------
      INTEGER i,ia,ib,ic,type
      REAL*8  furpar,qforce,alphar,twrtpi,gsrtal,xab,yab,zab,rsq,rsp
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,fsrtal(2),h,c1,c2
     &     ,c3,c4,qq,qq0,furpar0,corr0,corr,dcorr,en0,factl,factl0,facQ
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the intramolecular term ----------------------------------
*=======================================================================

      twrtpi=2.0d0/DSQRT(pi)
      fsrtal(1)=0.0D0
      fsrtal(2)=0.0D0
      en0=0.d0
      IF(erf_corr) THEN
#ifdef _OMP_
           if(omp_dynamic) call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC)  NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,list,ss_index,charge,lambdaq,lambdaq0,x0,y0,z0)
!$OMP& SHARED(intramol)
!$omp& SHARED(rlew,delew,erf_arr_corr,alphal,a1,a2,a3,a4,a5,qp,twrtpi)
!$OMP& REDUCTION(+:fsrtal,en0,fpx,fpy,fpz)  
#endif
        DO i=1,nlist
          ia=list(1,i)
          ib=list(2,i)
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
          
          type=ss_index(ia)
          facQ=list(3,i)
          if(intramol.and.list(3,i).ne.0) go to 101 ! skip all except on/off
          factl = (1.d0-dabs(lambdaq(ia)))*(1.-dabs(lambdaq(ib))) - facQ
          factl0 = (1.d0-dabs(lambdaq0(ia)))*(1.-dabs(lambdaq0(ib)))
     &         -facQ
          qq = factl*charge(ia)*charge(ib)
          qq0 = factl0*charge(ia)*charge(ib)
          xab=x0(ia)-x0(ib)
          yab=y0(ia)-y0(ib)
          zab=z0(ia)-z0(ib)
          rsq=xab*xab+yab*yab+zab*zab
          rsp=DSQRT(rsq)
          ic = 1+ (rsp-rlew)/delew  
          h = rsp-(rlew+(ic-1)*delew)
          c1=erf_arr_corr(1,ic) 
          c2=erf_arr_corr(2,ic) 
          c3=erf_arr_corr(3,ic) 
          c4=erf_arr_corr(4,ic) 
          corr =  (c1+h*(c2+h*(c3+h*c4/3)/2))
          corr0 = (c1+h*(c2+h*(c3+h*c4/3)/2))
          dcorr = (c2+h*(c3+h*c4/2) )/rsp
          alphar=alphal*rsp
          qt=1.0D0/(1.0e0+qp*alphar)
          expcst=exp(-alphar*alphar)
          erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x         expcst
          erfst=1.0D0-erfcst
          qforce= -qq*(dcorr +(erfst-twrtpi*alphar*expcst
     &         )/(rsp*rsq))
!$OMP ATOMIC
          fsrtal(type)=fsrtal(type)-qq*(erfst/rsp - corr)
          en0=en0 -qq0*(erfst/rsp - corr0)
          fpx(ia)=fpx(ia)+qforce*xab 
          fpy(ia)=fpy(ia)+qforce*yab
          fpz(ia)=fpz(ia)+qforce*zab
          fpx(ib)=fpx(ib)-qforce*xab
          fpy(ib)=fpy(ib)-qforce*yab
          fpz(ib)=fpz(ib)-qforce*zab
101       CONTINUE
        END DO
#ifdef _OMP_
!$OMP       END PARALLEL DO
#endif
      ELSE 
#ifdef _OMP_
           if(omp_dynamic)  call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO SCHEDULE(STATIC)  DEFAULT(PRIVATE) NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,list,ss_index,charge,lambdaq,lambdaq0,x0,y0,z0)
!$OMP& SHARED(alphal,a1,a2,a3,a4,a5,qp,twrtpi,intramol)
!$OMP& REDUCTION(+:fsrtal,en0,fpx,fpy,fpz)  
#endif
        DO i=1,nlist
          ia=list(1,i)
          ib=list(2,i)
          
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
          
          type=ss_index(ia)
          facQ=list(3,i)
          if(intramol.and.list(3,i).ne.0) go to 201
          factl = (1.d0-dabs(lambdaq(ia)))*(1.-dabs(lambdaq(ib))) -facQ
          factl0 = (1.d0-dabs(lambdaq0(ia)))*(1.-dabs(lambdaq0(ib)))
     &         -facQ
          furpar=factl*charge(ia)*charge(ib)
          furpar0=factl0*charge(ia)*charge(ib)
          xab=x0(ia)-x0(ib)
          yab=y0(ia)-y0(ib)
          zab=z0(ia)-z0(ib)
          rsq=xab*xab+yab*yab+zab*zab
          rsp=DSQRT(rsq)
          alphar=alphal*rsp
          qt=1.0D0/(1.0e0+qp*alphar)
          expcst=exp(-alphar*alphar)
          erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x         expcst
          erfst=1.0D0-erfcst
          qforce=-furpar*(erfst-twrtpi*alphar*expcst)/(rsp*rsq)
!$OMP ATOMIC
          fsrtal(type)=fsrtal(type)-furpar*erfst/rsp
          en0=en0-furpar0*erfst/rsp
          fpx(ia)=fpx(ia)+qforce*xab
          fpy(ia)=fpy(ia)+qforce*yab
          fpz(ia)=fpz(ia)+qforce*zab
          fpx(ib)=fpx(ib)-qforce*xab
          fpy(ib)=fpy(ib)-qforce*yab
          fpz(ib)=fpz(ib)-qforce*zab
201       CONTINUE
        END DO
#ifdef _OMP_
!$OMP       END PARALLEL DO
#endif
      END IF
#ifdef _OMP_ 
      if(omp_dynamic)  THEN 
        call omp_set_dynamic(.FALSE.)
        call OMP_SET_NUM_THREADS(nthr) !on exit restore thread numebr
      END IF
#endif
      fscnstr_slt=fsrtal(1)
      fscnstr_slv=fsrtal(2)
      dwrk=fsrtal(1)+fsrtal(2)-en0

*================= END OF EXECUTABLE STATEMENTS ========================
      RETURN
      END
