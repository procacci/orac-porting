      SUBROUTINE force(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,charge
     &     ,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &     ,ecc12,ecc6,alphal,erfc_spline,erfc_bin,erfc_arr
     &     ,mapnl,ngrp,grppt,uconf_slt,uconf_slv,uconf_ss,ucoul_slt
     &     ,ucoul_slv,ucoul_ss,fpx,fpy,fpz,stress,nnlpp,nnlpp1,rneigh
     &     ,rinn,rout,rtolinn,rtolout,gmass,iz,ncount1,rcut_corr,mass
     &     ,massinfty,ewald,lskip_ewald)


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

!****MultipleTimeScale Version*****P.Procacci-CECAM*********************
!
!!!   suboutine mts_forpp     (NON VECTOR version  )
!======================= DECLARATIONS ==================================
!****MultipleTimeScale Version*****P.Procacci-CECAM*********************

      use unit
      use parst
#ifdef PARALLEL
      use orac_mpi
#endif       

      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,ngrp,ma
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),groupp(*)
     &     ,atomp(*),atomg(*)
      INTEGER  nnlpp(*),nnlpp1(*)
      INTEGER mapnl(*),iz,ncount1
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*),ecc12(*),ecc6(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),fpx(*),fpy(*)
     &     ,fpz(*),gmass(*),erfc_arr(4,*),stress(3,3),mass(*),massinfty
      REAL*8  alphal,uconf_slt,uconf_slv,uconf_ss,ucoul_slt,ucoul_slv
     &     ,ucoul_ss,rneigh,rinn,rout,rtolinn,rtolout,erfc_bin,rcut_corr
      LOGICAL ewald,erfc_spline,lskip_ewald

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,jj,li,lj,lij,n,m,mm,noff,nbti,mapa,mapb,la,na,map,
     x     j1,i1,mp0,mp1,mbeg,ind,typei,typej,typeij,p1,p2,k,ic
      CHARACTER*80 errmsg
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x     xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,eps1
      REAL*8 ssvir,r6,r12,chrgei,auxa,auxb,c1,c2,c3,c4,h
      REAL*8 ucon,a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,conf,st(3,3)
      REAL*8 r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1
     &     ,arsout2,arsinn1,arsinn2,xmap0j,ymap0j,zmap0j,uconf(3)
     &     ,ucoul(3),st1,st2,st3,st4,st5,st6,st7,st8
     &     ,st9,xpcmp1,ypcmp1,zpcmp1,xxcm,yycm,zzcm,emvir,qfx,qfy,qfz
      REAL*8 xx,dxx,aux1,zero,one,two,three,rspi,derfcst,xd,yd,zd,xd1
     &     ,yd1,zd1,ucoula,uconfa,twoi,threei
      REAL*8 fpxt(nato),fpyt(nato),fpzt(nato)
#ifdef PARALLEL
      REAL *8 fpxtt(nato),fpytt(nato),fpztt(nato),ucoult(3),uconft(3)
     &     ,stresst(3,3)
      INTEGER ierr,root 
#endif
      INTEGER istride,itask,mp
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,      
     x     three=3.0D0,twoi=0.5D0,threei=1.0D0/3.0D0)

C------------------ DEFINITION OF A SCRATCH COMMON ---------------------

      INTEGER nb
      PARAMETER(nb=m1)
      INTEGER index(nb),indexa(nb),indexb(nb),index0(nb),
     &     bindex(nb),bindexa(nb),bindexb(nb),index1(nb),lcountb
      LOGICAL maplg(nb),mapag(nb)
      REAL*8 xmap0(tgroup),ymap0(tgroup),zmap0(tgroup)
      REAL*8 xmap3(tgroup),ymap3(tgroup),zmap3(tgroup),dr(tgroup)
      REAL*8 xmap1(tgroup),ymap1(tgroup),zmap1(tgroup),cmap2(tgroup)
     &     ,swrs(tgroup),dswrs(tgroup),erftbdns,xmap2(tgroup)
     &     ,ymap2(tgroup),zmap2(tgroup)
      REAL*8  rsqt(nb),rspt(nb),rsqit(nb),swrst(nb),ecc6t(nb),ecc12t(nb)
     &     ,charget(nb),conft(nb),coult(nb),ssvirt(nb)
      INTEGER typef(nb),j1t(nb),jt(nb),ifalse              

      COMMON /RAG1/ xmap3,ymap3,zmap3,xmap1,ymap1,zmap1,xmap2,ymap2
     &     ,zmap2,xmap0,ymap0,zmap0,cmap2,swrs,dswrs,bindex,bindexa
     &     ,index1,indexa,index,indexb,bindexb,lcountb,maplg,mapag

      REAL*8 PBC,x
      PBC(x)=DNINT(0.5D0*x)
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

C==================== EXECUTABLE STATEMENTS ============================
#ifdef PARALLEL
      itask=1+iproc_f
      istride=nproc_f
#else
      istride=1
      itask=1
#endif

      IF(nb.LT.nato) THEN
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
      n=0
      ncount1=0
      na=0


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

      DO j=1,nato
        maplg(j)=.TRUE.
      end do
#ifdef PARALLEL 
      if(iz.eq.1.and.rinn0.lt.rout) THEN ! only zeroes forces if not phony update call
        DO j=1,nato
          fpxt(i)=0.d0
          fpyt(i)=0.d0
          fpzt(i)=0.d0
          fpxtt(i)=0.d0
          fpytt(i)=0.d0
          fpztt(i)=0.d0
        end do
      endif
#endif
      DO j=1,ngrp
        mapag(j)=.TRUE.
      END DO

c==== start outer loop on groups
      twrtpi=2.0d0/DSQRT(pi)
      mbeg=1
      DO i=itask,ngrp,istride
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         m=nnlpp(1+n)
         map=0
         mp=0
         typei=ss_index(grppt(1,i))
c------------------------------------------------------------------------
c        build up neighbor list for the next inner shell and
c        set up group-group contact maps for:
c        interactions rinn0<r<rout (non switched) [INDEX0]
c             "       (rinn<r<rinn0.and.rout<r<rout0) (switched) [INDEX1]
c------------------------------------------------------------------------
         if(iz.eq.0) mbeg=i+1
         DO jj=mbeg,m
            j=nnlpp(jj+1+n)
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
               nnlpp1(1+map+ncount1)=j
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
2003        CONTINUE
         END DO
         IF(rinn.gt.1.d-5) then 
            nnlpp1(1+ncount1)=map
            ncount1=ncount1+map+1
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
               noff=mapnl(1+na)+1
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
            la=mapnl(1+na)
            p1=atomp(i1)
            xpcmp1=xpcm(p1)
            ypcmp1=ypcm(p1)
            zpcmp1=zpcm(p1)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.FALSE.
            END DO
            lcountb=0
c--         1-2 1-3 and 1-4 contacts are skipped
c--         switched group interactions are taken here...
            DO jj=1,mapb
               j=index1(jj)
               typej=ss_index(grppt(1,j))
               DO j1=grppt(1,j),grppt(2,j)
                 IF(maplg(j1).and.(mass(j1).lt.massinfty.or.mass(i1).lt
     &                .massinfty)) THEN
                     lcountb=lcountb+1
                     p2=atomp(j1)
                     xxcm = xpcmp1 - xpcm(p2) - xmap1(jj)
                     yycm = ypcmp1 - ypcm(p2) - ymap1(jj)
                     zzcm = zpcmp1 - zpcm(p2) - zmap1(jj)
                     xg=xpi-xp0(j1)-xmap1(jj)
                     yg=ypi-yp0(j1)-ymap1(jj)
                     zg=zpi-zp0(j1)-zmap1(jj)
                     xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                     yc=           co(2,2)*yg+co(2,3)*zg
                     zc=                      co(3,3)*zg
                     rsqt(lcountb)=xc*xc+yc*yc+zc*zc
                     rspt(lcountb)=DSQRT(rsqt(lcountb))
                     rsqit(lcountb)=1.0d0/rsqt(lcountb)
                     typef(lcountb)=typei+typej-1
                     swrst(lcountb)=swrs(jj)
                     ecc6t(lcountb)=ecc6(type(nbti,nbtype(j1)))
                     ecc12t(lcountb)=ecc12(type(nbti,nbtype(j1)))
                     j1t(lcountb)=jj
                     jt(lcountb)=j1
                     charget(lcountb)=charge(j1)
                  END IF
               END DO
            END DO
c----       compute switched forces:
c-----      calculate S(r)dV/dr and map sum_{k1k2} V_{k1,k2} onto cmap2
CDIR$    IVDEP
            DO jj=1,lcountb
               qforce =0.d0
               ifalse=0
               if(ecc6t(jj).gt.1.D-9) THEN ! Do LJ 
                 r6=rsqit(jj)*rsqit(jj)*rsqit(jj)
                 r12=r6*r6
                 ssvir=12.0d0*ecc12t(jj)*r12-6.0d0*ecc6t(jj)*r6
                 qforce=ssvir*rsqit(jj)*swrst(jj)
                 conft(jj)=ecc12t(jj)*r12-ecc6t(jj)*r6
                 ifalse=ifalse+1
               END IF
               IF(.not.lskip_ewald)THEN ! Do Ewald
                 rspqi=rsqit(jj)/rspt(jj)
                 alphar=alphal*rspt(jj)
                 qt=1.0d0/(1.0d0+qp*alphar)
                 expcst=DEXP(-alphar*alphar)
                 erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x                *qt+a1)*qt*expcst
                 furpar=chrgei*charget(jj)
                 coult(jj)=furpar*erfcst/rspt(jj)
                 qforce=qforce+furpar*(erfcst+twrtpi*alphar
     x              *expcst)*rspqi*swrst(jj)
                 ifalse=ifalse+1
               else if(.not.ewald) THEN  ! do standard coulomb
                 ssvirt(jj)=chrgei*charget(jj)/rsp
                 qforce=qforce+ssvirt(jj)*rsqi*swrst(jj)
                 ifalse=ifalse+1
               end if
               if(ifalse.ne.0) THEN 
                 emvir=qforce
                 fpx(i1)=fpx(i1)+qforce*xc
                 fpy(i1)=fpy(i1)+qforce*yc
                 fpz(i1)=fpz(i1)+qforce*zc
                 fpxt(jj)=-qforce*xc
                 fpyt(jj)=-qforce*yc
                 fpzt(jj)=-qforce*zc
                 st1 = st1+emvir*xc*xxcm
                 st2 = st2+emvir*xc*yycm
                 st3 = st3+emvir*xc*zzcm
                 st4 = st4+emvir*yc*xxcm
                 st5 = st5+emvir*yc*yycm
                 st6 = st6+emvir*yc*zzcm
                 st7 = st7+emvir*zc*xxcm
                 st8 = st8+emvir*zc*yycm
                 st9 = st9+emvir*zc*zzcm
               ENDIF
            END DO
            DO  jj=1,lcountb
              typeij=typef(jj)
              j1=j1t(jj)
              j=jt(jj)
#ifdef PARALLEL               
              fpxtt(j)=fpxtt(j)+fpxt(jj)
              fpytt(j)=fpytt(j)+fpyt(jj)
              fpztt(j)=fpztt(j)+fpzt(jj)
#else
              fpx(j)=fpx(j)+fpxt(jj)
              fpy(j)=fpy(j)+fpyt(jj)
              fpz(j)=fpz(j)+fpzt(jj)
#endif
              if(ecc6t(jj).gt.1.d-9) THEN 
                cmap2(j1)=cmap2(j1)+conft(jj)
                uconf(typeij)=uconf(typeij)+swrs(j1)*conft(jj)
              end if
              if(.not.lskip_ewald) THEN
                ucoul(typeij)=ucoul(typeij)+swrs(j1)*coult(jj)
                cmap2(j1)=cmap2(j1)+coult(jj)
              else if(.not.ewald) THEN 
                cmap2(j1)=cmap2(j1)+ssvir
                ucoul(typeij)=ucoul(typeij)+ssvir*swrst(jj)
              ENDIF
            END DO
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.TRUE.
            END DO
            noff=mapnl(1+na)+1
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
                  fpx(j1)=fpx(j1)-massj*xmap3(jj)
                  fpy(j1)=fpy(j1)-massj*ymap3(jj)
                  fpz(j1)=fpz(j1)-massj*zmap3(jj)
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
1002     noff=nnlpp(1+n)+1
c---     if IZ.EQ.1 (outer neighbor list on) do not rewind nnlww array
         if(iz.eq.1) n=n+noff
         IF(mpp.LT. ncount1) THEN
            errmsg=' IN MTS_FORPP : Dimensions of the NNLPP '/ /
     x           'are insufficient. ABORT !! '
            WRITE(kprint,'(a)') errmsg
            WRITE(kprint,'('' O1 ='',i7,'' NCOUNT = '',i7)')mpp,ncount1
            CALL xerror(errmsg,80,1,2)
         END IF
      END DO
      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)
      uconf_slt=uconf(1)
      uconf_ss =uconf(2)
      uconf_slv=uconf(3)
      IF(DABS(rout-rinn) .GT. 10D-5) THEN
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
      END IF
#ifdef PARALLEL
!     no forces are actually computed if force is a phony call for neighbor lists         
!     no further calculations are need; just exit 
      CALL MPI_BARRIER(newcomm,ierr) ! barrier at the end.
      if(iz.eq.0.or.rinn0.gt.rout) go to 101 
!     if beyond this point, forces have been computed and must be collected 
      root=0
      if(nproc_f.eq.1) THEN 
        write(kprint,*) "radii in force", rinn,rout,iproc_f,iproc_t
        do i=1,nato
          fpx(i)=fpx(i)+fpxt(i) 
          fpy(i)=fpy(i)+fpyt(i) 
          fpz(i)=fpz(i)+fpzt(i) 
          write(kprint,5000) i,fpx(i),fpxtt(i)
        end do
        goto 101                ! Only one task: no coomunications needed
      end if
      CALL MPI_BARRIER(newcomm,ierr) ! barrier at the end.
      write(kprint,*) "Forces before REDUCE",iproc_f,iproc_t
      do i=1,nato
        write(kprint,5000) i,fpx(i),fpxt(i)
      end do
!parallel sum of forces on J atoms (results on root of each communication group)
!message passing is here. Starting a message requires 1000 times the
!sending of data. The cost in MPI is starting a message and NOT passing
!the data. Hence the code below must minimize the CALL MPI_something 

      do i=1,nato 
        fpxtt(i)=0.0d0
        fpytt(i)=0.0d0
        fpztt(i)=0.0d0
      end do
      CALL MPI_REDUCE (fpxt,fpxtt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      CALL MPI_REDUCE (fpyt,fpytt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      CALL MPI_REDUCE (fpzt,fpztt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)

!now broadcast the results to all members of replica communication group.
      call MPI_BCAST(fpxtt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      call MPI_BCAST(fpytt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      call MPI_BCAST(fpztt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
! get the total force on i atoms (each proc doing its own job). proc i
! has the complete force on his own particles and zero otherwise. 
      write(kprint,*) "Forces AFTER REDUCE",iproc_f,iproc_t
      do i=itask,nato,istride
        fpx(i)=fpx(i)+fpxtt(i)
        fpy(i)=fpy(i)+fpytt(i)
        fpz(i)=fpz(i)+fpztt(i)
        write(kprint,5000) i,fpx(i),fpxtt(i)
5000    FORMAT(i10,2G15.5)
      end do
      
      if(nproc_f == 1) go to 101 !Only one task: no need for communications

!now sum all forces and send them to root in newcomm 
!root has now the complete force vector 
      do i=1,nato
        fpxtt(i)=0.d0
        fpytt(i)=0.d0
        fpztt(i)=0.d0
      end do
!!  NB: This iss done for debugging purposes: ONLY TWO CALL TO REDUCE AND !!! 
!!  BCAST ARE NEEDED
      CALL MPI_REDUCE (fpx,fpxtt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      CALL MPI_REDUCE (fpy,fpytt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      CALL MPI_REDUCE (fpz,fpztt,nato,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
    
 
    ! now broadcast new ftot from root to all within newcomm 
      call MPI_BCAST(fpxtt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      call MPI_BCAST(fpytt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      call MPI_BCAST(fpztt, nato, MPI_DOUBLE_PRECISION, 0, newcomm,
     &   ierr)
    
      do i=1,nato
        fpx(i)=fpxtt(i)
        fpy(i)=fpytt(i)
        fpz(i)=fpztt(i)
      end do

!  same stuff is done for energies and virial 
      CALL MPI_REDUCE (uconf,uconft,3,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      call MPI_BCAST(uconft, 3, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      do i=1,3
        uconf(i)=uconft(i)
      end do
      CALL MPI_REDUCE (ucoul,ucoult,3,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      call MPI_BCAST(ucoult, 3, MPI_DOUBLE_PRECISION, 0, newcomm,
     &     ierr)
      do i=1,3
        ucoul(i)=ucoult(i)
      end do
      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)

      uconf_slt=uconf(1)
      uconf_ss =uconf(2)
      uconf_slv=uconf(3)
      CALL MPI_REDUCE (stress,stresst,9,MPI_DOUBLE_PRECISION, MPI_SUM,
     &     root, newcomm,ierr)
      call MPI_BCAST(stresst, 9, MPI_DOUBLE_PRECISION, 0, newcomm, ierr)
      do i=1,3
        do j=1,3
          stress(i,j)=stresst(i,j)
        end do
      end do
101   CONTINUE
      write(kprint,*) "energies" 
      write(kprint,104) ucoul_slt,ucoul_ss,ucoul_slv
      write(kprint,104) uconf_slt,uconf_ss,uconf_slv
      write(kprint,*) "stress" 
      do i=1,3
        write(kprint,104) (stress(i,j),j=1,3)
104     format(3f12.5)
      end do
#else
      if(rinn0.lt.rout) THEN
      write(kprint,*) "Forces",rinn0,rout 
      do i=1,nato
        write(kprint,5000) i,fpx(i)
5000    FORMAT(i10,G15.5)
      end do
      write(kprint,*) "energies" 
      write(kprint,104) ucoul_slt,ucoul_ss,ucoul_slv
      write(kprint,104) uconf_slt,uconf_ss,uconf_slv
      write(kprint,*) "stress" 
      do i=1,3
        write(kprint,104) (stress(i,j),j=1,3)
104     format(3f12.5)
      end do
      endif
#endif
      RETURN

c========END OF 1 SUM FORCE ROUTINE================================

      END
