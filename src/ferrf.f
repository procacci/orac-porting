      SUBROUTINE ferrf(ss_index,alphal,charge,fudge,x0,y0,z0,iz,list
     &     ,nlist,fscnstr_slt,fscnstr_slv,fpx,fpy,fpz,erf_corr
     &     ,erf_arr_corr,delew,rlew,nbonds_added,mb)

************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald correction deriving from intramolecular        *
*     interaction.                                                     *
*                                                                      *
*---- Last update 06/12/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
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
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iz,nlist,list(2,*),ss_index(*),nbonds_added,mb
      REAL*8  charge(*),x0(*),y0(*),z0(*),fpx(mb),fpy(mb),fpz(mb)
     &     ,fscnstr_slt,fscnstr_slv,erf_arr_corr(4,*),delew,rlew 
      REAL*8  alphal,fudge
      LOGICAL erf_corr

*-------------------- LOCAL VARIABLES ----------------------------------
      INTEGER i,ia,ib,ic,type
      REAL*8  furpar,qforce,alphar,twrtpi,gsrtal,xab,yab,zab,rsq,rsp
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,fsrtal(2),h,c1,c2
     &     ,c3,c4,qq,corr,dcorr
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      SAVE fsrtal

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the intramolecular term ----------------------------------
*=======================================================================

      twrtpi=2.0d0/DSQRT(pi)
      fsrtal(1)=0.0
      fsrtal(2)=0.0
      IF(iz.EQ.1) THEN
         IF(erf_corr) THEN
#ifdef _OMP_
           if(omp_dynamic) call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC)  NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,nbonds_added,list,ss_index,charge,fudge,x0,y0,z0)
!$OMP& SHARED(rlew,delew,erf_arr_corr,alphal,a1,a2,a3,a4,a5,qp,twrtpi)
!$OMP& REDUCTION(+:fsrtal,fpx,fpy,fpz)  
#endif
            DO i=1,nlist-nbonds_added
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
               type=ss_index(ia)
               qq = charge(ia)*charge(ib)
               furpar=fudge*qq
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
               corr = qq*(c1+h*(c2+h*(c3+h*c4/3)/2))
               dcorr = qq*(c2+h*(c3+h*c4/2) )/rsp
               alphar=alphal*rsp
               qt=1.0D0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
               qforce= -dcorr - furpar*(erfst-twrtpi*alphar*expcst
     &              )/(rsp*rsq)
!$OMP ATOMIC
               fsrtal(type)=fsrtal(type)-(furpar*erfst/rsp - corr)
               fpx(ia)=fpx(ia)+qforce*xab 
               fpy(ia)=fpy(ia)+qforce*yab
               fpz(ia)=fpz(ia)+qforce*zab
               fpx(ib)=fpx(ib)-qforce*xab
               fpy(ib)=fpy(ib)-qforce*yab
               fpz(ib)=fpz(ib)-qforce*zab
            END DO
#ifdef _OMP_
!$OMP    END PARALLEL DO
#endif
         ELSE
#ifdef _OMP_
           if(omp_dynamic)  call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO SCHEDULE(STATIC)  DEFAULT(PRIVATE) NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,nbonds_added,list,ss_index,charge,fudge,x0,y0,z0)
!$OMP& SHARED(alphal,a1,a2,a3,a4,a5,qp,twrtpi)
!$OMP& REDUCTION(+:fsrtal,fpx,fpy,fpz)  
#endif
            DO i=1,nlist-nbonds_added
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=fudge*charge(ia)*charge(ib)
               xab=x0(ia)-x0(ib)
               yab=y0(ia)-y0(ib)
               zab=z0(ia)-z0(ib)
               rsq=xab*xab+yab*yab+zab*zab
               rsp=DSQRT(rsq)
               alphar=alphal*rsp
               qt=1.0D0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
               qforce=-furpar*(erfst-twrtpi*alphar*expcst)/(rsp*rsq)
!$OMP ATOMIC
               fsrtal(type)=fsrtal(type)-furpar*erfst/rsp
               fpx(ia)=fpx(ia)+qforce*xab
               fpy(ia)=fpy(ia)+qforce*yab
               fpz(ia)=fpz(ia)+qforce*zab
               fpx(ib)=fpx(ib)-qforce*xab
               fpy(ib)=fpy(ib)-qforce*yab
               fpz(ib)=fpz(ib)-qforce*zab
            END DO
#ifdef _OMP_
!$OMP       END PARALLEL DO
#endif
         END IF
      ELSE
         IF(erf_corr) THEN
#ifdef _OMP_
      if(omp_dynamic)  call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC)  NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,nbonds_added,list,ss_index,charge,x0,y0,z0)
!$OMP& SHARED(rlew,delew,erf_arr_corr,alphal,a1,a2,a3,a4,a5,qp,twrtpi)
!$OMP& REDUCTION(+:fsrtal)  
#endif
            DO i=1,nlist-nbonds_added
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=charge(ia)*charge(ib)
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
               corr = c1+h*(c2+h*(c3+h*c4/3)/2)
               alphar=alphal*rsp
               qt=1.0e0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
               erfst=1.0D0-erfcst
!$OMP ATOMIC
               fsrtal(type)=fsrtal(type)-furpar*(erfst/rsp - corr)
            END DO
#ifdef _OMP_
!$OMP       END PARALLEL DO
#endif
         ELSE
#ifdef _OMP_
           if(omp_dynamic) call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) NUM_THREADS(nthr1)
!$OMP& SHARED(nlist,nbonds_added,list,ss_index,charge,x0,y0,z0)
!$OMP& SHARED(alphal,a1,a2,a3,a4,a5,qp,twrtpi)
!$OMP& REDUCTION(+:fsrtal)  
#endif
            DO i=1,nlist-nbonds_added
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=charge(ia)*charge(ib)
               xab=x0(ia)-x0(ib)
               yab=y0(ia)-y0(ib)
               zab=z0(ia)-z0(ib)
               rsq=xab*xab+yab*yab+zab*zab
               rsp=DSQRT(rsq)
               alphar=alphal*rsp
               qt=1.0e0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
!$OMP ATOMIC
               fsrtal(type)=fsrtal(type)-furpar*erfst/
     x              rsp
            END DO
#ifdef _OMP_
!$OMP       END PARALLEL DO
#endif
         END IF
      END IF
#ifdef _OMP_ 
      if(omp_dynamic)  THEN 
        call omp_set_dynamic(.FALSE.)
        call OMP_SET_NUM_THREADS(nthr) !on exit restore thread numebr
      END IF
#endif
      fscnstr_slt=fsrtal(1)
      fscnstr_slv=fsrtal(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
