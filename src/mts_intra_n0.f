      subroutine mts_intra_n0(xp0,yp0,zp0,xcm,ycm,zcm,fpx,fpy,fpz
     &     ,ubond_slt,ubond_slv,ubend_slt,ubend_slv,uitors_slt
     &     ,uitors_slv,mb)

************************************************************************
*     FINTRAP stretching and bending routines for solute.
*     (used only with r-RESPA on)
************************************************************************

      use parst
      use cpropar
#ifdef _OMP_
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
      use unit
#endif       
      use omp_integr
      use rem
      
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER   mb
      REAL*8  xp0(*),yp0(*),zp0(*),xcm(*),ycm(*),zcm(*),ubend_slt
     &     ,ubond_slt,uitors_slt,ubend_slv,ubond_slv,uitors_slv
      REAL*8  fpx(mb),fpy(mb),fpz(mb)

*-------------------- LOCAL VARIABLES GROUP SCALING---------------------
#ifdef _OMP_
!     definitions below no longer needed: aux force; 
!     group scaling can't be done. 
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#else
      REAL*8  fppx(mb,2),fppy(mb,2),fppz(mb,2)
#endif
      integer itask
*-------------------- LOCAL VARIABLES ----------------------------------
      
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      REAL*8  xc,yc,zc,xd,yd,zd,sumx,sumy,sumz,st(3,3)
      REAL*8 seg_ene
      INTEGER i,j,i1,count,m,iat

      if(omp_timing) THEN 
        call timer(treal,tcpu1,tela1)
      end if

*=======================================================================
*----- Bonded interactions: stretching ------------------------------------
*=======================================================================

#ifdef _OMP_
      if(omp_dynamic)  call omp_set_dynamic(.TRUE.)
      ubond_slt=0.d0
      ubond_slv=0.d0
      ubend_slt=0.d0
      ubend_slv=0.d0
      uitors_slt=0.d0
      uitors_slv=0.d0
!$OMP    PARALLEL  PRIVATE(ITASK) NUM_THREADS(nthr1)
!$OMP&   REDUCTION(+:ubond_slt,ubond_slv,ubend_slt,ubend_slv)
!$OMP&   REDUCTION(+:uitors_slt,uitors_slv,seg_ene,fpx,fpy,fpz)
      itask=1+OMP_GET_THREAD_NUM()
#else
      itask=1
#endif
      IF(stretch) THEN
         CALL fpbond(ss_index,lstrtch,lstretch,xp0,yp0,zp0,potbo(1,2)
     &        ,potbo(1,1),ubond_slt,ubond_slv,fpx,fpy,fpz,0,mb,itask)
      END IF

*=======================================================================
*----- Bonded interactions: bending ------------------------------------
*=======================================================================

      if(bending) THEN 
        CALL fpbend(ss_index,lbndg,lbend,xp0,yp0,zp0
     &     ,potbe(1,2),potbe(1,1),potbe(1,4),potbe(1,3),potbe(1,5)
     &       ,potbe(1,6),potbe(1,7),ubend_slt,ubend_slv,fpx,fpy,fpz,mb
     &       ,itask)
      END IF

*=======================================================================
*----- Bonded interactions: improper torsions --------------------------
*=======================================================================

      if(rem_segment) then
        CALL fitors_rem(ss_index,litr,litor,xp0,yp0,zp0,potit(1,2),
     x       potit(1,1),potit(1,3),uitors_slt,uitors_slv,fpx,fpy,fpz
     &       ,remitor,mb,itask,seg_ene)
      else
        CALL fitors(ss_index,litr,litor,xp0,yp0,zp0,potit(1,2),
     x       potit(1,1),potit(1,3),uitors_slt,uitors_slv,fpx,fpy,fpz,mb
     &       ,itask)
      endif
*=======================================================================
*---- Compute stress tensor if coupling is by group --------------------
*=======================================================================

#ifdef _OMP_
!$OMP END PARALLEL
      if(omp_dynamic)  THEN 
        call omp_set_dynamic(.FALSE.)
        call OMP_SET_NUM_THREADS(nthr) !restore thread numebr
      END IF
#endif

      if(rem_segment) rem_n0_energy=seg_ene

!     write(6,'(2G15.5)') seg_ene,rem_n0_energy
c--   set to zero forces on fixed atoms
      do i=1,nfixed
        iat=ifixed(i)
        fpx(iat)=0.d0
        fpy(iat)=0.d0
        fpz(iat)=0.d0
      end do
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2) 
        time_n0=time_n0+tela2-tela1
        tcpu_n0=tcpu_n0+tcpu2-tcpu1
      end if
      RETURN
      END
