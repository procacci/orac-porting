c----------------------------------------------------------------------
      SUBROUTINE force(ss_index,xp0,yp0,zp0,xpg,ypg,zpg
     &     ,charge,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp
     &     ,co,alphal,mapnl,ngrp,grppt,uconf_slt,uconf_slv
     &     ,uconf_ss,ucoul_slt,ucoul_slv,ucoul_ss,fpx,fpy,fpz,stress
     &     ,nnlpp,nnlpp1,rneigh,rinn,rout,rtolinn,rtolout,gmass,iz
     &     ,ncount1,mass,massinfty,lskip_ewald,ewald,rmin,epsm)
c----------------------------------------------------------------------
c   Force routine (in vector shape) 
c---------------------------------------------------------------------- 
*****MultipleTimeScale Version*****P.Procacci-CECAM*********************
*                                                                      *
*     Compute the contributio from non-bonded interaction to          *
*     the forces and energies of the macromolecule. It does not        *
*     include code for hydrogen bonded interactions. Interactions      *
*     are switched with a group-group 3-spline function. Works both    *
*     for Ewald and NoEwald.                                           *
*                                                                      *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     CHARGE  :  List of atomic charges for the macromolecule.    (I)  *
*                >> real*8 CHARGE(NATO) <<                             *
*     NBTYPE  :  List of atomic types for the macromolecule.      (I)  *
*                >> real*8 NBTYPE(NATO) <<                             *
*     NATO    :  Number of atoms forming the macromolecule.       (I)  *
*     CO      :  Transformation matrix from box coordinates       (I)  *
*                to orthogonal frame.                                  *
*                >> real*8 CO(3,3) <<                                  *
*     SIGMA   :  List of L-J sigma parameters.                    (I)  *
*     EPSM    :  List of L-J epsilon parameters.                  (I)  *
*                >> real*8         <<                                  *
*     LKSIP_EWALD:  Logical parameter. If .TRUE. the electrostatic(I)  *
*                interaction (Ewald)  must be skipped                  *
*                >> logical*4 1 <<                                     *
*     EWALD   :  Logical parameter. If .TRUE. do Ewald elec.      (I)  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     MAPNL   :  Integer 1-2 and 1-3 list.                        (I)  *
*                >> integer   MAPNL(*) <<                              *
*     UCONF   :  Configurational LJ energy.                      (I/O) *
*     UCOUL   :  Coulombic energy.                               (I/O) *
*                                                                      *
*====================== WORK ARRAYS ===================================*
*                                                                      *
*     XMAP0   :  >> real*8  XMAP1(group)<<       Non-switched          *
*     YMAP0   :  >> real*8  YMAP1(group)<<       pbc vectors           *
*     ZMAP0   :  >> real*8  ZMAP1(group)<<                             *
*     XMAP1   :  >> real*8  XMAP1(group)<<       Switched pbc          *
*     YMAP1   :  >> real*8  YMAP1(group)<<       vectors               *
*     ZMAP1   :  >> real*8  ZMAP1(group)<<                             *
*     XMAP3   :  >> real*8  XMAP1(group)<<       switched              *
*     YMAP3   :  >> real*8  YMAP1(group)<<       distance vectors      *
*     ZMAP3   :  >> real*8  ZMAP1(group)<<                             *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************


C======================= DECLARATIONS ==================================

      use unit
      use parst
      
#ifdef _OMP_
      use omp_integr, only: nthr,m8t,mpp8,array_omp
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#endif
      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER iz,ma,nato,ngrp,groupp(*),atomp(*)
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),atomg(*)
#ifdef _OMP_
      INTEGER  nnlpp(mpp8,nthr),nnlpp1(mpp8,nthr),mapnl(m8t,nthr)
     &     ,ncount1(nthr)
      REAL*8 fpx(nato),fpy(nato),fpz(nato),fpxt(nato),fpyt(nato)
     &     ,fpzt(nato)
#else
      INTEGER nnlpp(*),nnlpp1(*)
      INTEGER mapnl(*),ncount1
      REAL*8 fpx(*),fpy(*),fpz(*)
#endif
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),gmass(*)
     &     ,stress(3,3),mass(*),massinfty
      REAL*8  alphal,uconf_slt,uconf_slv,uconf_ss,ucoul_slt,ucoul_slv
     &     ,ucoul_ss,rneigh,rinn,rout,rtolinn,rtolout
     &     ,rmin(*),epsm(*)
      LOGICAL lskip_ewald,ewald

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,jj,li,lj,lij,n,m,k,mm,noff,nbti,mapa,mapb,la,na,map,
     x     j1,i1,mp0,mp1,mbeg,mp,typei,typej,typeij,p1,p2,ic
      CHARACTER*80 errmsg
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x     xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,eps1
      REAL*8 ssvir,ssvirr,ssvirl,r6,r12,chrgei,aux1,auxa,auxb,drsq,c1,c2
     &     ,c3,c4,h
      REAL*8 ucon,a1,a2,a3,a4,a5,qp,qt,expcst,erfcst, st(3,3)
      REAL*8 rspqi,alphar,furpar,twrtpi,conf
      REAL*8 r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1
     &     ,arsout2,arsinn1,arsinn2,xmap0j,ymap0j,zmap0j,uconf(3)
     &     ,ucoul(3),st1,st2,st3,st4,st5,st6,st7,st8,st9,xpcmp1,ypcmp1
     &     ,zpcmp1,xxcm,yycm,zzcm,emvir,qfx,qfy,qfz,auxc
      INTEGER istride,itask,dimmax,ndims

!----------------WORK STORAGE AREA   540 K for all systems--------------
      REAL*8 xmap3(ngrpmax),ymap3(ngrpmax),zmap3(ngrpmax),dr(ngrpmax)
      REAL*8 xmap1(ngrpmax),ymap1(ngrpmax),zmap1(ngrpmax),cmap2(ngrpmax)
     &     ,xmap2(ngrpmax),ymap2(ngrpmax),zmap2(ngrpmax),swrs(ngrpmax)
     &     ,dswrs(ngrpmax)
      REAL*8  sigma,epsmin
      INTEGER bindex(nato),bindexa(nato),bindexb(nato),index1(nato)
     &     ,lcountb
      LOGICAL maplg(nato)
!---  Stoprage area for KNL architecture ----------------------------
      real*8     v_rmin(ngrpmax),v_epsm(ngrpmax),v_xpcm(ngrpmax)
     &     ,v_ypcm(ngrpmax),v_zpcm(ngrpmax),v_xmap1(ngrpmax)
     &     ,v_ymap1(ngrpmax),v_zmap1(ngrpmax),v_xp0(ngrpmax)
     &     ,v_yp0(ngrpmax),v_zp0(ngrpmax),v_swrs(ngrpmax)
     &     ,v_cmap2(ngrpmax),v_uconf(ngrpmax)
     &     ,v_ucoul(ngrpmax),v_charge(ngrpmax)
      logical ljj(ngrpmax)
#ifdef _OMP_
      REAL*8   v_fpxt(ngrpmax),v_fpyt(ngrpmax),v_fpzt(ngrpmax)
#else
      REAL*8   v_fpx(ngrpmax),v_fpy(ngrpmax),v_fpz(ngrpmax)
#endif      

!---------------- END OF WORK STORAGE AREA ----------------------------

      REAL*8 PBC,x
      PBC(x)=DNINT(0.5D0*x)
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer for all platforms but BGQ
#endif
#endif

C==================== EXECUTABLE STATEMENTS ============================

      IF(m1.LT.nato) THEN
        errmsg=' IN FNBOND : Dimensions of the work arrays insufficien
     xt. ABORTC '
        WRITE(6,'(a)') errmsg
        STOP
      END IF

*=======================================================================
*---- set up the big bunch of cut-off radii
*=======================================================================

c---  radius for the inner shell
      r2inn=rinn**2

c---  radius for the outer shell
      r2out=rout**2

c---  inner shell radius + healing lenght
      rinn0 = rinn+rtolinn
      r2inn0 =rinn0**2
      if(rinn.lt.1.d-5) rtolinn=0.d0

c---  outer shell radius +healing lenght
      rout0 = rout+rtolout
      r2out0 = rout0**2

c---- tolerance radius for the neighbor list of the inner shell
c---- to be used for shorter ranged interactions
      r2neigh=(rinn0+rneigh)**2

      twrtpi=2.0d0/DSQRT(pi)
      arsout1=rout0-3.0d0*rout
      arsout2=rtolout**3
      arsinn1=rinn0-3.0d0*rinn
      arsinn2=rtolinn**3
      st1=0.0D0
      st2=0.0D0
      st3=0.0D0
      st4=0.0D0
      st5=0.0D0
      st6=0.0D0
      st7=0.0D0
      st8=0.0D0
      st9=0.0D0
      DO i=1,3
        uconf(i)=0.0D0
        ucoul(i)=0.0D0
      END DO
#ifdef  _OMP_
      if(iz.eq.1.and.rinn0.lt.rout) THEN ! only zeroes forces if not phony update call
        DO j=1,nato
          fpxt(j)=0.d0
          fpyt(j)=0.d0
          fpzt(j)=0.d0
        end do
      endif
!$OMP  PARALLEL DEFAULT(PRIVATE) 
!$OMP& REDUCTION(+:fpxt,fpyt,fpzt)
!$OMP& REDUCTION(+:st1,st2,st3,st4,st5,st6,st7,st8,st9)
!$OMP& REDUCTION(+:ucoul,uconf) 
!$OMP& SHARED(rmin,epsm)
!$OMP& SHARED(xpg,ypg,zpg,ss_index,co,nnlpp,nnlpp1)
!$OMP& SHARED(grppt,xp0,yp0,zp0,nbtype,charge)
!$OMP& SHARED(atomg,fpx,fpy,fpz)
!$OMP& SHARED(a1,a2,a3,a4,a5,qp,r2inn,r2out,rinn0)
!$OMP& SHARED(r2inn0,rout0,r2out0,r2neigh)
!$OMP& SHARED(arsout1,arsout2,arsinn1,arsinn2,twrtpi,iz,gmass,mass)
!$OMP& SHARED(massinfty,atomp,xpcm,ypcm,zpcm)
!$OMP& SHARED(type,mapnl,nthr,m8t,mpp8)  
!$OMP& SHARED(nato,ma,ngrp,rinn,rout,ncount1,dimmax,ndims)
!$OMP& SHARED(lskip_ewald,ewald,alphal,istride,array_omp,ffield_type)
!$OMP& PRIVATE(n,na,itask,maplg,mbeg,xpgi,ypgi,zpgi,m)
!$OMP& PRIVATE(map,mp,typei,j,xgg,ygg,zgg,xmap0j,ymap0j,zmap0j)
!$OMP& PRIVATE(xc,yc,zc,drj,xmap1,ymap1,zmap1,xmap2,ymap2,zmap2)
!$OMP& PRIVATE(dr,index1,mapb,swrs,dswrs,rsp,auxa,auxb,cmap2)
!$OMP& PRIVATE(noff,xpi,ypi,zpi,nbti,chrgei,p1,la,mm)
!$OMP& PRIVATE(lcountb,xpcmp1,ypcmp1,zpcmp1,typej,bindex)
!$OMP& PRIVATE(bindexa,bindexb)
!$OMP& PRIVATE(sigma,epsmin,j1,typeij)
!$OMP& PRIVATE(lij,p2,xxcm,yycm,zzcm,xg,yg,zg,rsq,r6)
!$OMP& PRIVATE(ssvir,auxc,ucon)
!$OMP& PRIVATE(qforce,alphar,qt,expcst,erfcst,furpar)
!$OMP& PRIVATE(emvir)
!$OMP& PRIVATE(xmap3)
!$OMP& PRIVATE(ymap3)
!$OMP& PRIVATE(zmap3)
      n=0  ! why this must be here? 
      na=0
      itask=1+OMP_GET_THREAD_NUM()
!     write(6,*) "itask",itask,istride
      ncount1(itask)=0
      dimmax=mpp8
#else
      dimmax=mpp
      n=0
      ncount1=0
      na=0
#endif
      DO j=1,nato
        maplg(j)=.TRUE.
      END DO
c==== start outer loop on groups
      mbeg=1
      DO i=1,ngrp
#ifdef _OMP_
         if(array_omp(i,itask)) THEN
#endif
        xpgi=xpg(i)
        ypgi=ypg(i)
        zpgi=zpg(i)
#ifdef _OMP_
        m=nnlpp(1+n,itask)
#else
        m=nnlpp(1+n)
#endif
        map=0
        mp=0
        typei=ss_index(grppt(1,i))

c------------------------------------------------------------------------
c        build up neighbor list for the next inner shell and
c        set up group-group contact maps for:
c        interactions rinn0<r<rout (non switched) [INDEX0]
c             "       (rinn<r<rinn0.and.rout<r<rout0) (switched)
c             [INDEX1]
c------------------------------------------------------------------------
        if(iz.eq.0) mbeg=i+1
        DO jj=mbeg,m
#ifdef _OMP_
          j=nnlpp(jj+1+n,itask)
#else
          j=nnlpp(jj+1+n)
#endif
          xgg=xpgi-xpg(j)
          ygg=ypgi-ypg(j)
          zgg=zpgi-zpg(j)
          xmap0j=2.0D0*PBC(xgg)
          ymap0j=2.0D0*PBC(ygg)
          zmap0j=2.0D0*PBC(zgg)
          xgg=xgg-xmap0j
          ygg=ygg-ymap0j
          zgg=zgg-zmap0j
          xc=co(1,1)*xgg+co(1,2)*ygg+co(1,3)*zgg
          yc=            co(2,2)*ygg+co(2,3)*zgg
          zc=                        co(3,3)*zgg
          drj=xc**2+yc**2+zc**2
c---        do neighbor list for next inner shell
          if(r2inn.gt.1.0d-5.and.drj.LT.r2neigh) then
            map=map+1
#ifdef _OMP_
            nnlpp1(1+map+ncount1(itask),itask)=j
#else
            nnlpp1(1+map+ncount1)=j
#endif
          endif   
c---        skip mapping if outer shell is zero
          IF(rinn0.gt.rout) go to 2003
c---        map groups with "normal" + switched interactions
          IF(drj.GE.r2inn.AND.drj.LT.r2out0) THEN
            mp=mp+1
            xmap1(mp)=xmap0j
            ymap1(mp)=ymap0j
            zmap1(mp)=zmap0j
            xmap2(mp)=xc
            ymap2(mp)=yc
            zmap2(mp)=zc
            dr(mp)=drj
            index1(mp)=j
          ENDIF   
2003      CONTINUE
        END DO

        IF(rinn.gt.1.d-5) then 
#ifdef _OMP_
          nnlpp1(1+ncount1(itask),itask)=map
          ncount1(itask)=ncount1(itask)+map+1
#else
          nnlpp1(1+ncount1)=map
          ncount1=ncount1+map+1
#endif          
        endif

        mapb=mp
        
*-----   Computes switch. function and derivatives for all interacting 
*-----   Molecules and store them in swrs and dswrs
        DO jj=1,mapb
c----       if (rinn0 < r < rout ) the switching function is 1 
          IF(dr(jj).LT. r2out.and.dr(jj).GE.r2inn0) THEN
            swrs(jj)=1.0D0
            dswrs(jj)=0.0D0
          ELSE
c----          if (rout < r < rout0 ) the switching function is S(x) 
            IF(dr(jj).ge.r2out) then
              rsp=DSQRT(dr(jj))
              auxa=(arsout1+2.0d0*rsp)/arsout2
              auxb=rout0-rsp
              swrs(jj)=auxa*auxb**2
              dswrs(jj)=-2.0d0*auxa*auxb+2.0d0*auxb**2/arsout2
              dswrs(jj)=dswrs(jj)/rsp
c----          if (rinn < r < rinn0 ) the switching function is 1 -S(x) 
            ELSE
              rsp=DSQRT(dr(jj))
              auxa=(arsinn1+2.0d0*rsp)/arsinn2
              auxb=rinn0-rsp
              swrs(jj)=1.d0-auxa*auxb**2
              dswrs(jj)=2.0d0*auxa*auxb-2.0d0*auxb**2/arsinn2
              dswrs(jj)=dswrs(jj)/rsp
            ENDIF   
          END IF
        END DO

c---     set auxiliary array cmap2 for swich to zero
        DO jj=1,mapb
          cmap2(jj)=0.d0
        END DO
        IF(mapb.eq.0) THEN
          DO i1=grppt(1,i),grppt(2,i)
#ifdef _OMP_
            noff=mapnl(1+na,itask)+1
#else
            noff=mapnl(1+na)+1
#endif
            na=na+noff
          END DO
          go to 1002
        END IF
        DO i1=grppt(1,i),grppt(2,i)
          xpi=xp0(i1)
          ypi=yp0(i1)
          zpi=zp0(i1)
          nbti=nbtype(i1)
          chrgei=charge(i1)
          p1=atomp(i1)
          xpcmp1=xpcm(p1)
          ypcmp1=ypcm(p1)
          zpcmp1=zpcm(p1)
#ifdef _OMP_
          la=mapnl(1+na,itask)
          DO j=1,la
            mm=mapnl(j+1+na,itask)
            maplg(mm)=.FALSE.
          END DO
#else
          la=mapnl(1+na)
          DO j=1,la
            mm=mapnl(j+1+na)
            maplg(mm)=.FALSE.
          END DO
#endif          
          lcountb=0
c--         1-2 1-3 and 1-4 contacts are skipped
c--         switched group interactions are taken here...
          DO jj=1,mapb
            j=index1(jj)
            typej=ss_index(grppt(1,j))
            DO j1=grppt(1,j),grppt(2,j)
              IF(maplg(j1).and.(mass(j1).lt.massinfty.or.mass(i1).lt
     &             .massinfty)) THEN
                lcountb=lcountb+1
                bindex(lcountb)=j1
                bindexa(lcountb)=jj
                bindexb(lcountb)=typei+typej-1
                v_rmin(lcountb)=rmin(j1)
                v_epsm(lcountb)=epsm(j1)
                p2=atomp(j1)
                v_xpcm(lcountb)=xpcm(p2)                  
                v_ypcm(lcountb)=ypcm(p2)  
                v_zpcm(lcountb)=zpcm(p2)  
                v_xmap1(lcountb)=xmap1(jj)
                v_ymap1(lcountb)=ymap1(jj)
                v_zmap1(lcountb)=zmap1(jj)
                v_xp0(lcountb)=xp0(j1)
                v_yp0(lcountb)=yp0(j1)
                v_zp0(lcountb)=zp0(j1)
                v_charge(lcountb)=charge(j1)
                v_swrs(lcountb)=swrs(jj)
              END IF
            END DO
          END DO

c----       compute switched forces: 
c-----      calculate S(r)dV/dr and map sum_{k1k2} V_{k1,k2} onto cmap2
!CDIR$    IVDEP
          DO jj=1,lcountb
            ljj(jj)=.true.
            IF(ffield_type(1:4).eq."opls") THEN 
              sigma= 2.d0*(v_rmin(jj)*rmin(i1))**0.5
            ELSE
              sigma= (v_rmin(jj)+rmin(i1))
            end if
            epsmin = 4.d0*sqrt(epsm(i1)*v_epsm(jj))
            if(epsmin.lt.1.D-12.AND.lskip_ewald.and.ewald) THEN 
              ljj(jj)=.false.
              go   to 3001
            END IF
            xxcm = xpcmp1 - v_xpcm(jj) - v_xmap1(jj)
            yycm = ypcmp1 - v_ypcm(jj) - v_ymap1(jj)
            zzcm = zpcmp1 - v_zpcm(jj) - v_zmap1(jj)
            xg=xpi-v_xp0(jj)-v_xmap1(jj)
            yg=ypi-v_yp0(jj)-v_ymap1(jj)
            zg=zpi-v_zp0(jj)-v_zmap1(jj)
            xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
            yc=           co(2,2)*yg+co(2,3)*zg
            zc=                      co(3,3)*zg
            rsq=xc*xc+yc*yc+zc*zc
            auxa=1.d0/dsqrt(rsq)
            auxb = auxa*auxa*auxa
            if(epsmin.gt.1.D-12) THEN ! skip LJ is pot is ZERO
              r6=sigma**6.0*auxb**2
              r12=r6*r6
              ucon= epsmin*(r12-r6)
              v_cmap2(jj)=ucon
              v_uconf(jj)=v_swrs(jj)*ucon
              ssvir=epsmin*(12.d0*r12-6.0d0*r6)
              qforce=ssvir*auxa*auxa*v_swrs(jj)
              v_cmap2(jj)=ucon
              v_uconf(jj)=v_swrs(jj)*ucon
            ELSE
              qforce=0.d0
              auxc=0.d0
              v_uconf(jj)=0.d0
              v_cmap2(jj)=0.d0
            ENDIF
            if(.not.lskip_ewald) THEN !skip Ewald if pot is zero in r-shell 
              alphar=alphal/auxa
              qt=1.0d0/(1.0d0+qp*alphar)
              expcst=DEXP(-alphar*alphar)
              erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x             *qt+a1)*qt*expcst
              furpar=chrgei*v_charge(jj)
              v_ucoul(jj)=v_swrs(jj)*furpar*erfcst*auxa
              v_cmap2(jj)=v_cmap2(jj)+furpar*erfcst*auxa
              qforce=qforce+furpar*(erfcst+twrtpi*alphar
     x             *expcst)*auxb*v_swrs(jj)
            else if(.not.ewald) THEN  ! do bare coulomb 
              ssvir=chrgei*v_charge(jj)*auxa
              v_cmap2(jj)=v_cmap2(jj)+ssvir
              v_ucoul(jj)=ssvir*v_swrs(jj)
              qforce=qforce+ssvir*v_swrs(jj)*auxa*auxa
            else
              v_ucoul(jj)=0.0d0
            end if
            emvir=qforce
            fpx(i1)=fpx(i1)+qforce*xc
            fpy(i1)=fpy(i1)+qforce*yc
            fpz(i1)=fpz(i1)+qforce*zc
#ifdef _OMP_
            v_fpxt(jj)=-qforce*xc
            v_fpyt(jj)=-qforce*yc
            v_fpzt(jj)=-qforce*zc
#else
            v_fpx(jj)=-qforce*xc
            v_fpy(jj)=-qforce*yc
            v_fpz(jj)=-qforce*zc
#endif
            
            st1 = st1+emvir*xc*xxcm
            st2 = st2+emvir*xc*yycm
            st3 = st3+emvir*xc*zzcm
            st4 = st4+emvir*yc*xxcm
            st5 = st5+emvir*yc*yycm
            st6 = st6+emvir*yc*zzcm
            st7 = st7+emvir*zc*xxcm
            st8 = st8+emvir*zc*yycm
            st9 = st9+emvir*zc*zzcm
3001        CONTINUE
          END DO
! collect j forces and energies
          DO jj=1,lcountb
            if(ljj(jj)) THEN 
              j=bindex(jj)
              j1=bindexa(jj)
              typeij=bindexb(jj)
#ifdef _OMP_
              fpxt(j)=fpxt(j)+v_fpxt(jj)
              fpyt(j)=fpyt(j)+v_fpyt(jj)
              fpzt(j)=fpzt(j)+v_fpzt(jj)
#else
              fpx(j)=fpx(j)+v_fpx(jj)
              fpy(j)=fpy(j)+v_fpy(jj)
              fpz(j)=fpz(j)+v_fpz(jj)
#endif
              ucoul(typeij)=ucoul(typeij)+v_ucoul(jj)
              uconf(typeij)=uconf(typeij)+v_uconf(jj)
              cmap2(j1)=cmap2(j1)+v_cmap2(jj)
            end if
          END DO

#ifdef _OMP_
          la=mapnl(1+na,itask)
          DO j=1,la
            mm=mapnl(j+1+na,itask)
            maplg(mm)=.TRUE.
          END DO
          noff=mapnl(1+na,itask)+1
#else
          la=mapnl(1+na)
          DO j=1,la
            mm=mapnl(j+1+na)
            maplg(mm)=.TRUE.
          END DO
          noff=mapnl(1+na)+1
#endif
          na=na+noff
        END DO
        
c===     add the S'*V term to the atomic forces

        DO jj=1,mapb
          xmap3(jj)=-dswrs(jj)*cmap2(jj)*xmap2(jj)
          ymap3(jj)=-dswrs(jj)*cmap2(jj)*ymap2(jj)
          zmap3(jj)=-dswrs(jj)*cmap2(jj)*zmap2(jj)
        END DO

        if(mapb.ne.0) then
          DO i1=grppt(1,i),grppt(2,i)
            massi=gmass(i1)

            p1=atomp(i1)
            xpcmp1=xpcm(p1)
            ypcmp1=ypcm(p1)
            zpcmp1=zpcm(p1)

            do jj=1,mapb
              fpx(i1)=fpx(i1)+massi*xmap3(jj)
              fpy(i1)=fpy(i1)+massi*ymap3(jj)
              fpz(i1)=fpz(i1)+massi*zmap3(jj)

              st1=st1+massi*xmap3(jj)*xpcmp1
              st2=st2+massi*xmap3(jj)*ypcmp1
              st3=st3+massi*xmap3(jj)*zpcmp1
              st4=st4+massi*ymap3(jj)*xpcmp1
              st5=st5+massi*ymap3(jj)*ypcmp1
              st6=st6+massi*ymap3(jj)*zpcmp1
              st7=st7+massi*zmap3(jj)*xpcmp1
              st8=st8+massi*zmap3(jj)*ypcmp1
              st9=st9+massi*zmap3(jj)*zpcmp1

            end do
          end do
          do jj=1,mapb
            j=index1(jj)

            xxcm=xmap1(jj)
            yycm=ymap1(jj)
            zzcm=zmap1(jj)

            do j1=grppt(1,j),grppt(2,j)
              massj=gmass(j1)
#ifdef _OMP_
              fpxt(j1)=fpxt(j1)-massj*xmap3(jj)
              fpyt(j1)=fpyt(j1)-massj*ymap3(jj)
              fpzt(j1)=fpzt(j1)-massj*zmap3(jj)
#else
              fpx(j1)=fpx(j1)-massj*xmap3(jj)
              fpy(j1)=fpy(j1)-massj*ymap3(jj)
              fpz(j1)=fpz(j1)-massj*zmap3(jj)
#endif
              p2=atomp(j1)
              xpcmp1=-(xpcm(p2)+xxcm)*massj
              ypcmp1=-(ypcm(p2)+yycm)*massj
              zpcmp1=-(zpcm(p2)+zzcm)*massj
              st1=st1+xmap3(jj)*xpcmp1
              st2=st2+xmap3(jj)*ypcmp1
              st3=st3+xmap3(jj)*zpcmp1
              st4=st4+ymap3(jj)*xpcmp1
              st5=st5+ymap3(jj)*ypcmp1
              st6=st6+ymap3(jj)*zpcmp1
              st7=st7+zmap3(jj)*xpcmp1
              st8=st8+zmap3(jj)*ypcmp1
              st9=st9+zmap3(jj)*zpcmp1

            end do   
          end do
        end if
#ifdef _OMP_
1002    noff=nnlpp(1+n,itask)+1
        ndims=ncount1(itask)
#else
1002    noff=nnlpp(1+n)+1
        ndims=ncount1
#endif
c---     if IZ.EQ.1 (outer neighbor list on) do not rewind nnlww array
        if(iz.eq.1) n=n+noff
        IF(dimmax.LT.ndims) THEN
          errmsg=' IN MTS_FORPP : Dimensions of the NNLPP '/ /
     x         'are insufficient. ABORT !! '
          WRITE(kprint,'(a)') errmsg
          WRITE(kprint,'('' O1 ='',i7,'' NCOUNT = '',i7)') mpp,ncount1
          WRITE(kprint,'('' Arrived at atom = '',i7)') i
          CALL xerror(errmsg,80,1,2)
        END IF
#ifdef _OMP_
        ENDIF  !close task clause at the end of external loop
#endif
      END DO
#ifdef _OMP_
!$OMP END PARALLEL
      if(iz.eq.0.or.rinn0.ge.rout) goto 102 
      do i=1,nato
        fpx(i)=fpx(i)+fpxt(i)
        fpy(i)=fpy(i)+fpyt(i)
        fpz(i)=fpz(i)+fpzt(i)
      end do
#endif
      if(iz.eq.0.or.rinn0.ge.rout) goto 102 
      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)
      uconf_slt=uconf(1)
      uconf_ss =uconf(2)
      uconf_slv=uconf(3)
      st(1,1)=st1
      st(1,2)=st2
      st(1,3)=st3
      st(2,1)=st4
      st(2,2)=st5
      st(2,3)=st6
      st(3,1)=st7
      st(3,2)=st8
      st(3,3)=st9
      
      DO i=1,3
        DO j=1,3
          DO k=1,3
            stress(i,j) = stress(i,j)+st(i,k)*co(j,k)
          END DO
        END DO
      END DO

#ifdef _OMP_
      IF(nato.LT.0) THEN    ! to debug set nato.GT.0
        if(iz.eq.0.or.rinn0.ge.rout) goto 104 
        write(6,201) rinn
     &       ,rout
201     format(/50("="), /,10x, "R_in = ",f8.2,"    R_out = ",f8.2,/
     &       ,50("="),/9x,6x,"slt",12x,"slv",12x,"s-s")
        write(6,103) "coul",ucoul_slt,ucoul_slv,ucoul_ss
        write(6,103) "conf",uconf_slt,uconf_slv,uconf_ss
        DO i=1,3
          write(6,103) "stss",(stress(i,j),j=1,3)
        END DO
103     FORMAT(A4,5x,3G15.5)
        if(iz.eq.1) THEN 
          auxa=0.d0
          auxb=0.d0
          aux1=0.d0
          do i=1,nato
            auxa=auxa+fpx(i)**2
            auxb=auxb+fpy(i)**2
            aux1=aux1+fpz(i)**2
          end do
          write(6,77) i,auxa,auxb,aux1,xp0(i),yp0(i),zp0(i)
77        format(i10," ForcesI",3G15.5,/,10x,"Coord",3G15.5) 
          auxa=0.d0
          auxb=0.d0
          aux1=0.d0
          do i=1,nato
            auxa=auxa+fpxt(i)**2
            auxb=auxb+fpyt(i)**2
            aux1=aux1+fpzt(i)**2
          end do
          write(6,78) i,auxa,auxb,aux1,xp0(i),yp0(i),zp0(i)
78        format(i10," ForcesJ",3G15.5,/,10x,"Coord",3G15.5) 
        END IF
104     CONTINUE
      END iF
!     End of debug region-------------------------
!     End of debug region-------------------------
!     End of debug region-------------------------
#endif
102   CONTINUE
      RETURN
      END
