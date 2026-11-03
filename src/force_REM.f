      SUBROUTINE force_rem(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,charge
     &     ,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &     ,ecc12,ecc6,alphal,mapnl,ngrp,grppt,uconf_slt,uconf_slv
     &     ,uconf_ss,ucoul_slt,ucoul_slv,ucoul_ss,fpx,fpy,fpz,stress
     &     ,nnlpp,nnlpp1,rneigh,rinn,rout,rtolinn,rtolout,gmass,iz
     &     ,ncount1,mass,massinfty,seg_ene,rematom,lskip_ewald,ewald)

*****MultipleTimeScale Version*****P.Procacci-CECAM*********************
*                                                                      *
*     Compute the contribution from non-bonded interaction to          *
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
*     ECC12   :  List of L-J repulsive parameters.                (I)  *
*                >> real*8 ECC12(*) <<                                 *
*     ECC6    :  List of L-J attractive parameters.               (I)  *
*                >> real*8 ECC6(*) <<                                  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     MAPNL   :  Integer 1-2 and 1-3 list.                        (I)  *
*                >> integer   MAPNL(*) <<                              *
*     UCONF   :  Configurational energy.                         (I/O) *
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
      use rem
      use cpropar, only:intramol,ljintra
      
#ifdef _OMP_
      use omp_integr, only: nthr,m8t,mpp8,array_omp
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#endif
      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,ngrp,ma
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),groupp(*)
     &     ,atomp(*),atomg(*)
#ifdef _OMP_
      INTEGER  nnlpp(mpp8,nthr),nnlpp1(mpp8,nthr),mapnl(m8t,nthr)
     &     ,ncount1(nthr)
      REAL*8 fpx(nato),fpy(nato),fpz(nato),fpxt(nato),fpyt(nato)
     &     ,fpzt(nato)
#else
      INTEGER nnlpp(*),nnlpp1(*),mapnl(*),ncount1
      REAL*8 fpx(*),fpy(*),fpz(*)
#endif
      INTEGER  iz,ntap
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*),ecc12(*),ecc6(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),gmass(*)
     &     ,stress(3,3),mass(*),massinfty
      REAL*8  alphal,uconf_slt,uconf_slv,uconf_ss,ucoul_slt,ucoul_slv
     &     ,ucoul_ss,rneigh,rinn,rout,rtolinn,rtolout
!!
      REAL*8 seg_ene
      LOGICAL lskip_ewald,ewald
!! 
      LOGICAL rematom(*)

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,jj,li,lj,lij,n,m,mm,noff,nbti,mapa,mapb,la,na,map,
     x     j1,i1,mp0,mp1,mbeg,ind,typei,typej,typeij,p1,p2,k,ic
      CHARACTER*80 errmsg
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x     xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,eps1
      REAL*8 ssvir,r6,r12,chrgei,auxa,auxb,c1,c2,c3,c4,h
      REAL*8 ucon,a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,conf,coul,st(3,3)
      REAL*8 r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1
     &     ,arsout2,arsinn1,arsinn2,xmap0j,ymap0j,zmap0j,uconf(3)
     &     ,ucoul(3),st1,st2,st3,st4,st5,st6,st7,st8,st9,xpcmp1,ypcmp1
     &     ,zpcmp1,xxcm,yycm,zzcm,emvir,qfx,qfy,qfz
      REAL*8 xx,dxx,aux1,zero,one,two,three,rspi,derfcst,xd,yd,zd,xd1
     &     ,yd1,zd1,ucoula,uconfa,twoi,threei
      INTEGER istride,itask,dimmax,ndims
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,      
     x     three=3.0D0,twoi=0.5D0,threei=1.0D0/3.0D0)
C------------------ WORK STORAGE AREA -----------------------------
      INTEGER index0(ngrpmax),index1(ngrpmax),lcountb
      LOGICAL maplg(nato),mapag(nato)
      REAL*8 xmap0(ngrpmax),ymap0(ngrpmax),zmap0(ngrpmax)
      REAL*8 xmap3(ngrpmax),ymap3(ngrpmax),zmap3(ngrpmax)
      REAL*8 xmap1(ngrpmax),ymap1(ngrpmax),zmap1(ngrpmax),cmap2(ngrpmax)
     &     ,swrs(ngrpmax),dswrs(ngrpmax),erftbdns,xmap2(ngrpmax)
     &     ,ymap2(ngrpmax),zmap2(ngrpmax),charge_old(nato)

      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif
      REAL*8 PBC,x
      PBC(x)=DNINT(0.5D0*x)


C==================== EXECUTABLE STATEMENTS ============================
      IF(m1.LT.nato) THEN
        errmsg=' IN FNBOND : Dimensions of the work arrays insufficien
     xt. ABORTC '
        WRITE(6,'(a)') errmsg
        STOP
      END IF
      
!     if (lskip_force) the REM "segment" has modified intranmol
!     nonbonded potential:
!     a gas-phase REM will run  with LJ only
!     (ljintra=.true.) or with both charge and lj turned off 

      IF(LSKIP_FORCE) THEN
        if(.not.ljintra)  then
          RETURN
        ELSE
          do j=1,nato
            charge_old(j)=charge(j)  ! save charges 
            charge(j)=0.d0      ! set charges for intra segment to zero 
          END DO
        endif
      ENDIF
!!

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
      seg_ene = 0.d0      
#ifdef  _OMP_
      if(iz.eq.1.and.rinn0.lt.rout) THEN ! only zeroes forces if not phony update call
        DO j=1,nato
          fpxt(j)=0.d0
          fpyt(j)=0.d0
          fpzt(j)=0.d0
        end do
      endif
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& REDUCTION(+:fpxt,fpyt,fpzt)
!$OMP& REDUCTION(+:st1,st2,st3,st4,st5,st6,st7,st8,st9)
!$OMP& REDUCTION(+:ucoul,uconf,seg_ene) 
!$OMP& SHARED(xpg,ypg,zpg,ss_index,co,nnlpp,nnlpp1)
!$OMP& SHARED(grppt,xp0,yp0,zp0,nbtype,charge)
!$OMP& SHARED(atomg,fpx,fpy,fpz)
!$OMP& SHARED(a1,a2,a3,a4,a5,qp,r2inn,r2out,rinn0)
!$OMP& SHARED(r2inn0,rout0,r2out0,r2neigh)
!$OMP& SHARED(arsout1,arsout2,arsinn1,arsinn2,twrtpi,iz,gmass,mass)
!$OMP& SHARED(massinfty,ecc6,ecc12,atomp,xpcm,ypcm,zpcm)
!$OMP& SHARED(type,mapnl,nthr,m8t,mpp8,array_omp)  
!$OMP& SHARED(nato,ma,ngrp,rinn,rout,ncount1,dimmax,ndims)
!$OMP& SHARED(lskip_ewald,ewald,alphal,istride)
!$OMP& SHARED(rematom,rem_factor,rem_segkind,kprint)
!$OMP& PRIVATE(n,na,itask,maplg,mapag,mbeg,xpgi,ypgi,zpgi,m)
!$OMP& PRIVATE(map,mp0,mp1,typei,xgg,ygg,zgg,xmap0j,ymap0j,zmap0j) 
!$OMP& PRIVATE(xc,yc,zc,drj,xmap0,ymap0,zmap0,index0,xmap1,ymap1,zmap1)
!$OMP& PRIVATE(xmap2,ymap2,zmap2,index1,rsp,auxa,auxb,swrs,dswrs,mapa)
!$OMP& PRIVATE(mapb,cmap2,noff,xpi,ypi,zpi,nbti,chrgei)
!$OMP& PRIVATE(la,mm,p1,xpcmp1,ypcmp1,zpcmp1,j1,typej,typeij,xd,yd,zd)
!$OMP& PRIVATE(xd1,yd1,zd1,uconfa,ucoula,lij,p2,xxcm,yycm,zzcm,xg,yg,zg)
!$OMP& PRIVATE(rsq,rsqi,qforce,r6,r12,ssvir,conf,rspqi,alphar,qt,expcst)
!$OMP& PRIVATE(erfcst,furpar,aux1,emvir,qfx,qfy,qfz,xmap3,ymap3,zmap3)
!$OMP& PRIVATE(massi,massj,errmsg,coul)
      n=0  
      na=0
      itask=1+OMP_GET_THREAD_NUM()
!      write(6,*) "itask",itask,istride
      ncount1(itask)=0
      dimmax=mpp8
#else
      dimmax=mpp
      n=0
      na=0
      ncount1=0
#endif
      DO j=1,nato
        maplg(j)=.TRUE.
      END DO
      DO j=1,ngrp
        mapag(j)=.TRUE.
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
        mp0=0
        mp1=0
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
          IF(rinn0.gt.rout.or.iz.eq.0) go to 2003
c---        map groups with "normal" (non switched) interactions
          IF (drj.GE.r2inn0.AND.drj.LT.r2out) THEN
            mp0=mp0+1
            xmap0(mp0)=xmap0j
            ymap0(mp0)=ymap0j
            zmap0(mp0)=zmap0j
            index0(mp0)=j
c---        map groups in outer and inner "switched" shell
          ELSE
c----          if (rout < r < rout0 ) the switching function is S(x) 
            IF (drj.GE.r2out.and.drj.lt.r2out0) THEN
              mp1=mp1+1
              xmap1(mp1)=xmap0j
              ymap1(mp1)=ymap0j
              zmap1(mp1)=zmap0j
              xmap2(mp1)=xc
              ymap2(mp1)=yc
              zmap2(mp1)=zc

              index1(mp1)=j
              rsp=DSQRT(drj)
              auxa=(arsout1+2.0d0*rsp)/arsout2
              auxb=rout0-rsp
              swrs(mp1)=auxa*auxb**2
              dswrs(mp1)=-2.0d0*auxa*auxb+2.0d0*auxb**2/arsout2
              dswrs(mp1)=dswrs(mp1)/rsp
            ELSE IF(drj.GE.r2inn.and.drj.LT.r2inn0) then
c----          if (rinn < r < rinn0 ) the switching function is 1 - S(x)
              mp1=mp1+1
              xmap1(mp1)=xmap0j
              ymap1(mp1)=ymap0j
              zmap1(mp1)=zmap0j
              xmap2(mp1)=xc
              ymap2(mp1)=yc
              zmap2(mp1)=zc

              index1(mp1)=j
              rsp=DSQRT(drj)
              auxa=(arsinn1+2.0d0*rsp)/arsinn2
              auxb=rinn0-rsp
              swrs(mp1)=1.d0-auxa*auxb**2
              dswrs(mp1)=2.0d0*auxa*auxb-2.0d0*auxb**2/arsinn2
              dswrs(mp1)=dswrs(mp1)/rsp
            ENDIF   
          ENDIF
2003      continue
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

c---    this is the number of normal (non switched) group interactions
        mapa=mp0

c---    this is the number of switched [S or (1-S)] group interactions
        mapb=mp1
        
c---    set auxiliary array cmap2 for swich to zero
        DO jj=1,mapb
          cmap2(jj)=0.d0
        END DO
        IF(mapa.eq.0.and.mapb.eq.0) THEN
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
#ifdef _OMP_
          la=mapnl(1+na,itask)
          DO j=1,la
            mm=mapnl(j+1+na,itask)
            maplg(mm)=.FALSE.
            mapag(atomg(mm))=.FALSE.
          END DO
#else
          la=mapnl(1+na)
          DO j=1,la
            mm=mapnl(j+1+na)
            maplg(mm)=.FALSE.
            mapag(atomg(mm))=.FALSE.
          END DO
#endif          
          p1=atomp(i1)
          xpcmp1=xpcm(p1)
          ypcmp1=ypcm(p1)
          zpcmp1=zpcm(p1)

!----       compute forces
!         4 inner loup can be alternatively open: 
!         check whether first atom of the group is non bonded; if so
!         proceed in inner atomic j-loop; else open the inner j-loop but
!         check if every j atom is non-bonded to i. For each of 
!         thse possibilities, two separate contact group-maps - for
!         switched
!         and non switched interactions - are unfolded.  
          

!         This is the group map for non switched interactions 
          DO jj=1,mapa
            j1=index0(jj)
            typej=ss_index(grppt(1,j1))
            typeij=typei+typej-1
            xd=-xmap0(jj)
            yd=-ymap0(jj)
            zd=-zmap0(jj)
            xd1=xpi+xd
            yd1=ypi+yd
            zd1=zpi+zd

            xd=xpcmp1+xd
            yd=ypcmp1+yd
            zd=zpcmp1+zd

            uconfa=0.0D0
            ucoula=0.0D0
            IF(mapag(j1)) THEN
              DO j=grppt(1,j1),grppt(2,j1)
                lij=type(nbti,nbtype(j))
                if(lskip_ewald.and.(ecc6(lij).LT.1D-12).and.
     &                 ewald) goto 3001
                p2=atomp(j)
                xxcm = xd - xpcm(p2)
                yycm = yd - ypcm(p2)
                zzcm = zd - zpcm(p2)
                xg=xd1-xp0(j)
                yg=yd1-yp0(j)
                zg=zd1-zp0(j)
                xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                yc=           co(2,2)*yg+co(2,3)*zg
                zc=                      co(3,3)*zg
                rsq=xc*xc+yc*yc+zc*zc
                rsqi=1.0d0/rsq
                qforce=0.d0
                conf=0.d0
                if(ecc6(lij).GT.1D-12) THEN
                  r6=rsqi*rsqi*rsqi
                  r12=r6*r6
                  ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)
     &                 *r6
                  qforce=ssvir*rsqi
                  conf=ecc12(lij)*r12-ecc6(lij)*r6
                  uconfa=uconfa+conf
                end if
                if(.not.lskip_ewald) THEN
                  rsp=DSQRT(rsq)
                  rspqi=rsqi/rsp
                  alphar=alphal*rsp
                  qt=1.0d0/(1.0d0+qp*alphar)
                  expcst=DEXP(-alphar*alphar)
                  erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x                 *qt+a1)*qt*expcst
                  furpar=chrgei*charge(j)
                  coul=furpar*erfcst/rsp 
                  ucoula=ucoula+coul
                  aux1  = furpar*(erfcst+twrtpi*alphar*expcst)
     &                 *rspqi
                  qforce=qforce+aux1
                else if (.not.ewald) THEN 
                  rsp=DSQRT(rsq)
                  coul=chrgei*charge(j)/rsp
                  ucoul(typeij)=ucoul(typeij)+coul
                  qforce=qforce+coul*rsqi
                end if
                if(seg_select(rem_segkind,rematom,i1,j)) then
                  seg_ene = seg_ene + conf + coul
                  qforce = qforce * rem_factor(3)
                endif
                emvir = qforce
                fpx(i1)=fpx(i1)+qforce*xc
                fpy(i1)=fpy(i1)+qforce*yc
                fpz(i1)=fpz(i1)+qforce*zc
#ifdef _OMP_
                fpxt(j)=fpxt(j)-qforce*xc
                fpyt(j)=fpyt(j)-qforce*yc
                fpzt(j)=fpzt(j)-qforce*zc
#else
                fpx(j)=fpx(j)-qforce*xc
                fpy(j)=fpy(j)-qforce*yc
                fpz(j)=fpz(j)-qforce*zc
#endif
                qfx=emvir*xc
                qfy=emvir*yc
                qfz=emvir*zc
                st1 = st1+qfx*xxcm
                st2 = st2+qfx*yycm
                st3 = st3+qfx*zzcm
                st4 = st4+qfy*xxcm
                st5 = st5+qfy*yycm
                st6 = st6+qfy*zzcm
                st7 = st7+qfz*xxcm
                st8 = st8+qfz*yycm
                st9 = st9+qfz*zzcm
3001            CONTINUE
              END DO
            ELSE
              DO j=grppt(1,j1),grppt(2,j1)
                IF(maplg(j).and.(mass(j).lt.massinfty.or
     &               .mass(i1).lt.massinfty)) THEN
                  lij=type(nbti,nbtype(j))
                  if(lskip_ewald.and.(ecc6(lij).LT.1D-12).and.  
     &                 ewald) goto 3002
                  p2=atomp(j)
                  xxcm = xd - xpcm(p2)
                  yycm = yd - ypcm(p2)
                  zzcm = zd - zpcm(p2)
                  xg=xd1-xp0(j)
                  yg=yd1-yp0(j)
                  zg=zd1-zp0(j)
                  xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                  yc=           co(2,2)*yg+co(2,3)*zg
                  zc=                      co(3,3)*zg
                  rsq=xc*xc+yc*yc+zc*zc
                  rsqi=1.0d0/rsq
                  qforce=0.d0
                  conf=0.d0
                  if(ecc6(lij).GT.1D-12) THEN
                    r6=rsqi*rsqi*rsqi
                    r12=r6*r6
                    ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij
     &                   )*r6
                    qforce=ssvir*rsqi
                    conf=ecc12(lij)*r12-ecc6(lij)*r6
                    uconfa=uconfa+conf
                  end if
                  if(.not.lskip_ewald) THEN 
                    rsp=DSQRT(rsq)
                    rspqi=rsqi/rsp
                    alphar=alphal*rsp
                    qt=1.0d0/(1.0d0+qp*alphar)
                    expcst=DEXP(-alphar*alphar)
                    erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x                   *qt+a1)*qt*expcst
                    furpar=chrgei*charge(j)
                    coul=furpar*erfcst/rsp
                    ucoula=ucoula+coul
                    aux1  = furpar*(erfcst+twrtpi*alphar
     &                   *expcst)*rspqi
                    qforce=qforce+aux1
                  else if (.not.ewald) THEN 
                    rsp=DSQRT(rsq)
                    coul=chrgei*charge(j)/rsp
                    ucoul(typeij)=ucoul(typeij)+coul
                    qforce=qforce+coul*rsqi
                  end if
                  if(seg_select(rem_segkind,rematom,i1,j))
     &                 then
                    seg_ene = seg_ene + conf + coul
                    qforce = qforce * rem_factor(3)
                  endif
                  emvir=qforce
                  fpx(i1)=fpx(i1)+qforce*xc
                  fpy(i1)=fpy(i1)+qforce*yc
                  fpz(i1)=fpz(i1)+qforce*zc
#ifdef _OMP_
                  fpxt(j)=fpxt(j)-qforce*xc
                  fpyt(j)=fpyt(j)-qforce*yc
                  fpzt(j)=fpzt(j)-qforce*zc
#else
                  fpx(j)=fpx(j)-qforce*xc
                  fpy(j)=fpy(j)-qforce*yc
                  fpz(j)=fpz(j)-qforce*zc
#endif
                  qfx=emvir*xc
                  qfy=emvir*yc
                  qfz=emvir*zc
                  st1 = st1+qfx*xxcm
                  st2 = st2+qfx*yycm
                  st3 = st3+qfx*zzcm
                  st4 = st4+qfy*xxcm
                  st5 = st5+qfy*yycm
                  st6 = st6+qfy*zzcm
                  st7 = st7+qfz*xxcm
                  st8 = st8+qfz*yycm
                  st9 = st9+qfz*zzcm
3002            CONTINUE
                END IF
              END DO
            END IF
            ucoul(typeij)=ucoul(typeij)+ucoula
            uconf(typeij)=uconf(typeij)+uconfa
          END DO
c----     compute switched forces:
c-----    calculate S(r)dV/dr and map sum_{k1k2} V_{k1,k2} onto cmap2

          DO jj=1,mapb
            j1=index1(jj)
            typej=ss_index(grppt(1,j1))
            typeij=typei+typej-1
            IF(mapag(j1)) THEN
              DO j=grppt(1,j1),grppt(2,j1)
                lij=type(nbti,nbtype(j))
                if(lskip_ewald.and.(ecc6(lij).LT.1D-12).and.
     &                 ewald) goto 3003
                p2=atomp(j)
                xxcm = xpcmp1 - xpcm(p2) - xmap1(jj)
                yycm = ypcmp1 - ypcm(p2) - ymap1(jj)
                zzcm = zpcmp1 - zpcm(p2) - zmap1(jj)
                xg=xpi-xp0(j)-xmap1(jj)
                yg=ypi-yp0(j)-ymap1(jj)
                zg=zpi-zp0(j)-zmap1(jj)
                xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                yc=           co(2,2)*yg+co(2,3)*zg
                zc=                      co(3,3)*zg
                rsq=xc*xc+yc*yc+zc*zc
                rsqi=1.0d0/rsq
                qforce=0.d0
                conf=0.d0
                if(ecc6(lij).GT.1D-12) THEN
                  r6=rsqi*rsqi*rsqi
                  r12=r6*r6
                  ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                  qforce=ssvir*rsqi*swrs(jj)
                  conf=ecc12(lij)*r12-ecc6(lij)*r6
                  if(seg_select(rem_segkind,rematom,i1,j)) then
                    cmap2(jj)=cmap2(jj) + conf*rem_factor(3)
                  else
                    cmap2(jj)=cmap2(jj) + conf
                  endif
                  uconf(typeij)=uconf(typeij)+swrs(jj)*conf
                end if
                if(.not.lskip_ewald) THEN
                  rsp=DSQRT(rsq)
                  rspqi=rsqi/rsp
                  alphar=alphal*rsp
                  qt=1.0d0/(1.0d0+qp*alphar)
                  expcst=DEXP(-alphar*alphar)
                  erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
                  furpar=chrgei*charge(j)
                  coul=furpar*erfcst/rsp
                  ucoul(typeij)=ucoul(typeij) + coul*swrs(jj)
                  if(seg_select(rem_segkind,rematom,i1,j)) then
                    cmap2(jj)=cmap2(jj)+coul*rem_factor(3)
                  else
                    cmap2(jj)=cmap2(jj) + coul
                  endif
                  aux1=furpar*(erfcst+twrtpi*alphar*expcst)
     &                 *rspqi*swrs(jj)
                  qforce=qforce+aux1
                else if (.not.ewald) THEN
                  coul=chrgei*charge(j)/rsp
                  ucoul(typeij)=ucoul(typeij)+coul*swrs(jj)
                  qforce=qforce+coul*rsqi*swrs(jj)
                  if(seg_select(rem_segkind,rematom,i1,j)) then
                    cmap2(jj)=cmap2(jj) + coul *rem_factor(3)
                  else
                    cmap2(jj)=cmap2(jj) + coul
                  endif
                end if
                if(seg_select(rem_segkind,rematom,i1,j)) then
                  seg_ene = seg_ene + swrs(jj)*conf+swrs(jj)*coul
                  qforce = qforce * rem_factor(3)
                endif
                emvir=qforce
                fpx(i1)=fpx(i1)+qforce*xc
                fpy(i1)=fpy(i1)+qforce*yc
                fpz(i1)=fpz(i1)+qforce*zc
#ifdef _OMP_
                fpxt(j)=fpxt(j)-qforce*xc
                fpyt(j)=fpyt(j)-qforce*yc
                fpzt(j)=fpzt(j)-qforce*zc
#else
                fpx(j)=fpx(j)-qforce*xc
                fpy(j)=fpy(j)-qforce*yc
                fpz(j)=fpz(j)-qforce*zc
#endif
                qfx=emvir*xc
                qfy=emvir*yc
                qfz=emvir*zc
                st1 = st1+qfx*xxcm
                st2 = st2+qfx*yycm
                st3 = st3+qfx*zzcm
                st4 = st4+qfy*xxcm
                st5 = st5+qfy*yycm
                st6 = st6+qfy*zzcm
                st7 = st7+qfz*xxcm
                st8 = st8+qfz*yycm
                st9 = st9+qfz*zzcm
3003          CONTINUE
              END DO
            ELSE
              DO j=grppt(1,j1),grppt(2,j1)
                IF(maplg(j).and.(mass(j).lt.massinfty.or
     &               .mass(i1).lt.massinfty)) THEN
                  lij=type(nbti,nbtype(j))
                  if(lskip_ewald.and.(ecc6(lij).LT.1D-12).and.
     &                 ewald) goto 3004
                  p2=atomp(j)
                  xxcm = xpcmp1 - xpcm(p2) - xmap1(jj)
                  yycm = ypcmp1 - ypcm(p2) - ymap1(jj)
                  zzcm = zpcmp1 - zpcm(p2) - zmap1(jj)
                  xg=xpi-xp0(j)-xmap1(jj)
                  yg=ypi-yp0(j)-ymap1(jj)
                  zg=zpi-zp0(j)-zmap1(jj)
                  xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                  yc=           co(2,2)*yg+co(2,3)*zg
                  zc=                      co(3,3)*zg
                  rsq=xc*xc+yc*yc+zc*zc
                  rsqi=1.0d0/rsq
                  qforce=0.d0
                  conf=0.d0
                  if(ecc6(lij).GT.1D-12) THEN
                    r6=rsqi*rsqi*rsqi
                    r12=r6*r6
                    ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                    qforce=ssvir*rsqi*swrs(jj)
                    conf=ecc12(lij)*r12-ecc6(lij)*r6
                    if(seg_select(rem_segkind,rematom,i1,j))then
                      cmap2(jj)=cmap2(jj) + conf *rem_factor(3)
                    else
                      cmap2(jj)=cmap2(jj) + conf
                    endif
                    uconf(typeij)=uconf(typeij)+swrs(jj)*conf
                  end if
                  if(.not.lskip_ewald) THEN
                    rsp=DSQRT(rsq)
                    rspqi=rsqi/rsp
                    alphar=alphal*rsp
                    qt=1.0d0/(1.0d0+qp*alphar)
                    expcst=DEXP(-alphar*alphar)
                    erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x                   *qt+a1)*qt*expcst
                    furpar=chrgei*charge(j)
                    coul=furpar*erfcst/rsp
                    ucoul(typeij)=ucoul(typeij)+swrs(jj)*coul
                    if(seg_select(rem_segkind,rematom,i1,j)) then
                      cmap2(jj)=cmap2(jj) + coul*rem_factor(3)
                    else
                      cmap2(jj)=cmap2(jj) + coul
                    endif
                    aux1=furpar*(erfcst+twrtpi*alphar*expcst)
     &                   *rspqi*swrs(jj)
                    qforce=qforce+aux1
                  else if (.not.ewald) THEN
                      coul=chrgei*charge(j)/rsp
                      ucoul(typeij)=ucoul(typeij)+coul*swrs(jj)
                      qforce=qforce+coul*rsqi*swrs(jj)
                      if(seg_select(rem_segkind,rematom,i1,j)) then
                        cmap2(jj)=cmap2(jj) + coul *rem_factor(3)
                      else
                        cmap2(jj)=cmap2(jj) + coul
                      endif
                    end if
                  if(seg_select(rem_segkind,rematom,i1,j))
     &                 then
                    seg_ene = seg_ene + swrs(jj)*conf+swrs(jj)*coul
                    qforce = qforce * rem_factor(3)
                  endif
                  emvir=qforce
                  fpx(i1)=fpx(i1)+qforce*xc
                  fpy(i1)=fpy(i1)+qforce*yc
                  fpz(i1)=fpz(i1)+qforce*zc
#ifdef _OMP_
                  fpxt(j)=fpxt(j)-qforce*xc
                  fpyt(j)=fpyt(j)-qforce*yc
                  fpzt(j)=fpzt(j)-qforce*zc
#else
                  fpx(j)=fpx(j)-qforce*xc
                  fpy(j)=fpy(j)-qforce*yc
                  fpz(j)=fpz(j)-qforce*zc
#endif
                  qfx=emvir*xc
                  qfy=emvir*yc
                  qfz=emvir*zc
                  st1 = st1+qfx*xxcm
                  st2 = st2+qfx*yycm
                  st3 = st3+qfx*zzcm
                  st4 = st4+qfy*xxcm
                  st5 = st5+qfy*yycm
                  st6 = st6+qfy*zzcm
                  st7 = st7+qfz*xxcm
                  st8 = st8+qfz*yycm
                  st9 = st9+qfz*zzcm
3004              CONTINUE
                END IF
              END DO
            END IF 
          END DO
#ifdef _OMP_
          la=mapnl(1+na,itask)

          DO j=1,la
            mm=mapnl(j+1+na,itask)
            maplg(mm)=.TRUE.
            mapag(atomg(mm))=.TRUE.
          END DO
          noff=mapnl(1+na,itask)+1
#else
          la=mapnl(1+na)
          DO j=1,la
            mm=mapnl(j+1+na)
            maplg(mm)=.TRUE.
            mapag(atomg(mm))=.TRUE.
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
        IF(dimmax.LT. ndims) THEN
          errmsg=' IN MTS_FORPP : Dimensions of the NNLPP '/ /
     x         'are insufficient. ABORT !! '
          WRITE(kprint,'(a)') errmsg
          WRITE(kprint,'('' O1 ='',i7,'' NCOUNT = '',i7)')mpp,ncount1
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
!     Debug writing if (nato.lt.0) do nothing-------------------------
!     Debug writing if (nato.lt.0) do nothing-------------------------
!     Debug writing if (nato.lt.0) do nothing-------------------------
#ifdef _OMP_
      IF(nato.LT.0) THEN
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
      if(ljintra.and.lskip_force) THEN
        do j=1,nato
          charge(j)=charge_old(j)
        END DO
      ENDIF
102   CONTINUE
      RETURN
      END
