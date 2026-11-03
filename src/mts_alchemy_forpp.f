c----------------------------------------------------------------------
      SUBROUTINE mts_alchemy_forpp(ss_index,xp0,yp0,zp0,xpg,ypg,zpg
     &     ,charge,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp
     &     ,co,ecc12,ecc6,ewald,alphal,erfc_spline,erfc_bin,erfc_arr
     &     ,mapnl,ngrp,grppt,uconf_slt,uconf_slv,uconf_ss,ucoul_slt
     &     ,ucoul_slv,ucoul_ss,fpx,fpy,fpz,stress,nnlpp,nnlpp1,rneigh
     &     ,rinn,rout,rtolinn,rtolout,gmass,iz,ncount1,rcut_corr,mass
     &     ,massinfty,lambda,lambda0,lambdaq,lambdaq0,rmin,epsm
     &     ,dwrk)
c----------------------------------------------------------------------
c   Force routine (in vector shape) for alchemical transformation 
c   with soft-core GROMOS potential (Procacci 2013)  
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
*     ECC12   :  List of L-J repulsive parameters.                (I)  *
*                >> real*8 ECC12(*) <<                                 *
*     ECC6    :  List of L-J attractive parameters.               (I)  *
*                >> real*8 ECC6(*) <<                                  *
*     EWALD   :  Logical parameter. If .TRUE. the electrostatic   (I)  *
*                interaction is compute with Ewald.                    *
*                >> logical*4 EWALD <<                                 *
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
      
      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER ma,nato,ngrp,groupp(*),atomp(*)
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),atomg(*)
      INTEGER nnlpp(*),nnlpp1(*)
      INTEGER mapnl(*),iz,ncount1
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*),ecc12(*),ecc6(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),fpx(*),fpy(*)
     &     ,fpz(*),gmass(*),erfc_arr(4,*),stress(3,3),mass(*),massinfty
      REAL*8  alphal,uconf_slt,uconf_slv,uconf_ss,ucoul_slt,ucoul_slv
     &     ,ucoul_ss,rneigh,rinn,rout,rtolinn,rtolout,erfc_bin,rcut_corr
     &     ,lambda(*),lambda0(*),lambdaq(*),lambdaq0(*),rmin(*),epsm(*)
     &     ,dwrk
      LOGICAL ewald,erfc_spline

      INTEGER nb
      PARAMETER(nb=m1)
      LOGICAL maplg(nb)

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,jj,li,lj,lij,n,m,k,mm,noff,nbti,mapa,mapb,la,na,map,
     x        j1,i1,mp0,mp1,mbeg,mp,typei,typej,typeij,p1,p2,ic
      CHARACTER*80 errmsg
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x       xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,eps1
      REAL*8 ssvir,ssvirr,ssvirl,r6,r12,chrgei,auxa,auxb,drsq,c1,c2,c3
     &     ,c4,h
      REAL*8 ucon,a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,swrs(tgroup)
     &     ,dswrs(tgroup),st(3,3)
      REAL*8 rspqi,alphar,furpar,twrtpi,conf,furpar0,auxc0
      REAL*8 r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1
     &     ,arsout2,arsinn1,arsinn2,xmap0j,ymap0j,zmap0j,uconf(3)
     &     ,ucoul(3),st1,st2,st3,st4,st5,st6,st7,st8,st9,xpcmp1,ypcmp1
     &     ,zpcmp1,xxcm,yycm,zzcm,emvir,qfx,qfy,qfz,alch_fact(nb,2)
     &     ,alch_factq(nb,2),alch_fact0(nb,2),alch_factq0(nb,2)
     &     ,auxc

c----------------- SCRATCH COMMON BLOCK --------------------------------

      REAL*8 xmap3(tgroup),ymap3(tgroup),zmap3(tgroup),dr(tgroup)
      REAL*8 xmap1(tgroup),ymap1(tgroup),zmap1(tgroup),cmap2(tgroup)
     &     ,xmap2(tgroup),ymap2(tgroup),zmap2(tgroup),erftbdns
      REAL*8  gp,gp1,gp2,gr,gl,sigma,epsmin,apnd
      INTEGER bindex(nb),bindexa(nb),bindexb(nb),index1(nb),lcountb

      REAL*8  gp0,gp10,gp20,furpar1,auxc1,uconf0(3),ucoul0(3)
     &     ,ucon0,ssvir0
      COMMON /RAG1/ xmap3,ymap3,zmap3,xmap1,ymap1,zmap1,cmap2,dr,xmap2
     &     ,ymap2,zmap2,bindex,bindexa,bindexb,index1,lcountb,maplg


      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp,apnd/0.3275911d0,0.5d0/
      INCLUDE 'pbc.h'

C==================== EXECUTABLE STATEMENTS ============================

      IF(nb.LT.nato) THEN
          errmsg=' IN FNBOND : Dimensions of the work arrays insufficien
     xt. ABORTC '
          WRITE(6,'(a)') errmsg
          STOP
      END IF

*=======================================================================
*---- set up the big bunch of cut-off radii
*=======================================================================

      IF(erfc_spline) erftbdns=1.0D0/erfc_bin
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
         uconf0(i)=0.0D0
         ucoul0(i)=0.0D0
      END DO

c========mts_forpp_EWALD=================================================
c---       take EWALD SUM 
c========================================================================

      DO j=1,nato
         maplg(j)=.TRUE.
      END DO
      if(.not.ewald) go to 2005
c==== start outer loop on groups
      twrtpi=2.0d0/DSQRT(pi)
      mbeg=1
      DO i=1,ngrp
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
                     bindex(lcountb)=j1
                     bindexa(lcountb)=jj
                     bindexb(lcountb)=typei+typej-1
                  END IF
               END DO
            END DO
c--         computes the old and new alchemical factors (non vector loop) 
            DO jj=1,lcountb
              j=bindex(jj)
! LENNARD-JONES shifted  potential between normal and alchemical (on/off) atoms   
              if(lambda(i1).eq.0.d0.and.lambda(j).ne.0.d0) THEN 
                alch_fact(jj,1)=dabs(lambda(j))
                alch_fact(jj,2)=1.d0-dabs(lambda(j))
                alch_fact0(jj,1)=dabs(lambda0(j))
                alch_fact0(jj,2)=1.d0-dabs(lambda0(j))
              ELSE if(lambda(i1).ne.0.d0.and.lambda(j).eq.0.d0) THEN  
                alch_fact(jj,1)=dabs(lambda(i1))
                alch_fact(jj,2)=1.d0-dabs(lambda(i1))
                alch_fact0(jj,1)=dabs(lambda0(i1))
                alch_fact0(jj,2)=1.d0-dabs(lambda0(i1))
!             normal potential if lambda's are both equal to zero. 
              ELSE if (lambda(i1).eq.0.d0.and.lambda(j).eq.0.d0) THEN   
                alch_fact(jj,1)=0.d0 
                alch_fact(jj,2)=1.d0
                alch_fact0(jj,1)=0.d0 
                alch_fact0(jj,2)=1.d0
!             alchemical on atoms do not feel alchemical off atoms
              ELSE if ((lambda(i1).ne.0.d0.and.lambda(j).ne.0.d0).  
     &         AND.(lambda(i1).ne.lambda(j))) THEN     
                alch_fact(jj,1)=0.d0 
                alch_fact(jj,2)=0.d0
                alch_fact0(jj,1)=0.d0 
                alch_fact0(jj,2)=0.d0
!             on atoms (or off atoms) interact via the normal potential
              ELSE if ((lambda(i1).ne.0.d0.and.lambda(j).ne.0.d0).  
     &          AND.(lambda(i1).eq.lambda(j))) THEN   
                alch_fact(jj,1)=0.d0 
                alch_fact(jj,2)=1.d0
                alch_fact0(jj,1)=0.d0 
                alch_fact0(jj,2)=1.d0
              END IF
!CHARGES      shifted  potential between normal and alchemical (on/off) atoms   
              if(lambdaq(i1).eq.0.d0.and.lambdaq(j).ne.0.d0) THEN 
                alch_factq(jj,1)=dabs(lambdaq(j))
                alch_factq(jj,2)=1.d0-dabs(lambdaq(j))
                alch_factq0(jj,1)=dabs(lambdaq0(j))
                alch_factq0(jj,2)=1.d0-dabs(lambdaq0(j))
              ELSE if(lambdaq(i1).ne.0.d0.and.lambdaq(j).eq.0.d0) THEN  
                alch_factq(jj,1)=dabs(lambdaq(i1))
                alch_factq(jj,2)=1.d0-dabs(lambdaq(i1))
                alch_factq0(jj,1)=dabs(lambdaq0(i1))
                alch_factq0(jj,2)=1.d0-dabs(lambdaq0(i1))
!             normal potential if lambda's are zero 
              ELSE if (lambdaq(i1).eq.0.d0.and.lambdaq(j).eq.0.d0) THEN   
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=1.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=1.d0
!             alchemical on atoms do not feel alchemical off atoms
              ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(j).ne.0.d0).  
     &         AND.(lambdaq(i1).ne.lambdaq(j))) THEN     
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=0.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=0.d0
!             on atoms (or off atoms) interact via the normal potential
              ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(j).ne.0.d0).  
     &          AND.(lambdaq(i1).eq.lambdaq(j))) THEN   
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=1.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=1.d0
              END IF
            END DO

c----       compute switched forces: 
c-----      calculate S(r)dV/dr and map sum_{k1k2} V_{k1,k2} onto cmap2
CDIR$    IVDEP
            DO jj=1,lcountb
               j=bindex(jj)
               sigma= (rmin(j)+rmin(i1))
               epsmin = 4.d0*sqrt(epsm(i1)*epsm(j))
               j1=bindexa(jj)
               typeij=bindexb(jj)
               lij=type(nbti,nbtype(j))

               p2=atomp(j)
               xxcm = xpcmp1 - xpcm(p2) - xmap1(j1)
               yycm = ypcmp1 - ypcm(p2) - ymap1(j1)
               zzcm = zpcmp1 - zpcm(p2) - zmap1(j1)

               xg=xpi-xp0(j)-xmap1(j1)
               yg=ypi-yp0(j)-ymap1(j1)
               zg=zpi-zp0(j)-zmap1(j1)
               xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
               yc=           co(2,2)*yg+co(2,3)*zg
               zc=                      co(3,3)*zg
               rsq=xc*xc+yc*yc+zc*zc
               auxa=dsqrt(rsq)
               r6=rsq*rsq*rsq/sigma**6.
               gp=(apnd*alch_fact(jj,1) + r6)
               gp1=1.d0/gp
               gp2=gp1*gp1 
               gp0=(apnd*alch_fact0(jj,1) + r6)
               gp10=1.d0/gp0
               gp20=gp10*gp10 
               gr=6.d0*r6/auxa
               ssvir=alch_fact(jj,2)*epsmin*(2.d0*gp1*gp2 - gp2) 
               ssvirr=gr*ssvir
               auxc=epsmin*(gp2 - gp1)
               auxc0=epsmin*(gp20 - gp10)
               ucon=alch_fact(jj,2)*auxc
               ucon0=alch_fact0(jj,2)*auxc0
               auxb = 1./(auxa*auxa*auxa)
               qforce=ssvirr*swrs(j1)/auxa
               cmap2(j1)=cmap2(j1)+ucon
               uconf(typeij)=uconf(typeij)+swrs(j1)*ucon
               uconf0(typeij)=uconf0(typeij)+swrs(j1)*ucon0
               alphar=alphal*auxa
               qt=1.0d0/(1.0d0+qp*alphar)
               expcst=DEXP(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x              *qt+a1)*qt*expcst
               furpar=alch_factq(jj,2)*chrgei*charge(j)
               furpar0=alch_factq0(jj,2)*chrgei*charge(j)
               ucoul(typeij)=ucoul(typeij)+swrs(j1)*furpar*erfcst/auxa
               ucoul0(typeij)=ucoul0(typeij)+swrs(j1)*furpar0*erfcst
     &              /auxa
               cmap2(j1)=cmap2(j1)+furpar*erfcst/auxa
               qforce=qforce+furpar*(erfcst+twrtpi*alphar
     x              *expcst)*auxb*swrs(j1)
               emvir=qforce
               fpx(i1)=fpx(i1)+qforce*xc
               fpy(i1)=fpy(i1)+qforce*yc
               fpz(i1)=fpz(i1)+qforce*zc
               fpx(j)=fpx(j)-qforce*xc
               fpy(j)=fpy(j)-qforce*yc
               fpz(j)=fpz(j)-qforce*zc

               st1 = st1+emvir*xc*xxcm
               st2 = st2+emvir*xc*yycm
               st3 = st3+emvir*xc*zzcm
               st4 = st4+emvir*yc*xxcm
               st5 = st5+emvir*yc*yycm
               st6 = st6+emvir*yc*zzcm
               st7 = st7+emvir*zc*xxcm
               st8 = st8+emvir*zc*yycm
               st9 = st9+emvir*zc*zzcm

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
      dwrk=0.d0 
c--   compute the work in the alchemical transofrmation
      do i=1,3
        auxa=ucoul(i)-ucoul0(i)
        auxb=uconf(i)-uconf0(i)
        dwrk=dwrk+auxa+auxb
      end do
        
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

      RETURN

c========END OF EWALD SUM FORCE ROUTINE================================

2005  CONTINUE

c=======mts_forpp NOEWALD================================================
c---       do NOT take EWALD SUM 
c========================================================================

      mbeg=1
      DO i=1,ngrp
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
            IF(rinn0.gt.rout) go to 3003
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
3003        CONTINUE
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
            go to 2002
         END IF
         DO i1=grppt(1,i),grppt(2,i)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)

            p1=atomp(i1)
            xpcmp1=xpcm(p1)
            ypcmp1=ypcm(p1)
            zpcmp1=zpcm(p1)

            nbti=nbtype(i1)
            chrgei=charge(i1)

            la=mapnl(1+na)
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
                     bindex(lcountb)=j1
                     bindexa(lcountb)=jj
                     bindexb(lcountb)=typei+typej-1
                  END IF
               END DO
            END DO

c--         computes the alchemical factors (non vector loop) 
            DO jj=1,lcountb
              j=bindex(jj)
!             LENNARD JONES
              if(lambda(i1).eq.0.d0.and.lambda(j).ne.0.d0) THEN 
                alch_fact(jj,1)=dabs(lambda(j))
                alch_fact(jj,2)=1.d0-dabs(lambda(j))
                alch_fact0(jj,1)=dabs(lambda0(j))
                alch_fact0(jj,2)=1.d0-dabs(lambda0(j))
              ELSE if(lambda(i1).ne.0.d0.and.lambda(j).eq.0.d0) THEN  
                alch_fact(jj,1)=dabs(lambda(i1))
                alch_fact(jj,2)=1.d0-dabs(lambda(i1))
                alch_fact0(jj,1)=dabs(lambda0(i1))
                alch_fact0(jj,2)=1.d0-dabs(lambda0(i1))
              ELSE if (lambda(i1).eq.0.d0.and.lambda(j).eq.0.d0) THEN   
                alch_fact(jj,1)=0.d0 
                alch_fact(jj,2)=1.d0
                alch_fact0(jj,1)=0.d0 
                alch_fact0(jj,2)=1.d0
              ELSE if ((lambda(i1).ne.0.d0.and.lambda(j).ne.0.d0).  
     &         AND.(lambda(i1).ne.lambda(j))) THEN     
                alch_fact(jj,1)=0.d0 
                alch_fact(jj,2)=0.d0
                alch_fact0(jj,1)=0.d0 
                alch_fact0(jj,2)=0.d0
               ELSE if ((lambda(i1).ne.0.d0.and.lambda(j).ne.0.d0).  
     &         AND.(lambda(i1).eq.lambda(j))) THEN   
                alch_fact(jj,1)=0.d0
                alch_fact(jj,2)=1.d0
                alch_fact0(jj,1)=0.d0
                alch_fact0(jj,2)=1.d0
              END IF
!             CHARGES
              if(lambdaq(i1).eq.0.d0.and.lambdaq(j).ne.0.d0) THEN 
                alch_factq(jj,1)=dabs(lambdaq(j))
                alch_factq(jj,2)=1.d0-dabs(lambdaq(j))
                alch_factq0(jj,1)=dabs(lambdaq0(j))
                alch_factq0(jj,2)=1.d0-dabs(lambdaq0(j))
              ELSE if(lambdaq(i1).ne.0.d0.and.lambdaq(j).eq.0.d0) THEN  
                alch_factq(jj,1)=dabs(lambdaq(i1))
                alch_factq(jj,2)=1.d0-dabs(lambdaq(i1))
                alch_factq0(jj,1)=dabs(lambdaq0(i1))
                alch_factq0(jj,2)=1.d0-dabs(lambdaq0(i1))
!             normal potential if lambda's are zero 
              ELSE if (lambdaq(i1).eq.0.d0.and.lambdaq(j).eq.0.d0) THEN   
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=1.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=1.d0
!             alchemical on atoms do not feel alchemical off atoms
              ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(j).ne.0.d0).  
     &         AND.(lambdaq(i1).ne.lambdaq(j))) THEN     
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=0.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=0.d0
!             on atoms (or off atoms) interact via the normal potential
              ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(j).ne.0.d0).  
     &          AND.(lambdaq(i1).eq.lambdaq(j))) THEN   
                alch_factq(jj,1)=0.d0 
                alch_factq(jj,2)=1.d0
                alch_factq0(jj,1)=0.d0 
                alch_factq0(jj,2)=1.d0
              END IF
            END DO
c----       compute switched forces
CDIR$    IVDEP
            DO jj=1,lcountb
               j=bindex(jj)
               sigma= (rmin(j)+rmin(i1))
               epsmin = 4.d0*sqrt(epsm(i1)*epsm(j))
               j1=bindexa(jj)
               typeij=bindexb(jj)
               lij=type(nbti,nbtype(j))
               
               p2=atomp(j)
               xxcm = xpcmp1 - xpcm(p2) - xmap1(j1)
               yycm = ypcmp1 - ypcm(p2) - ymap1(j1)
               zzcm = zpcmp1 - zpcm(p2) - zmap1(j1)

               xg=xpi-xp0(j)-xmap1(j1)
               yg=ypi-yp0(j)-ymap1(j1)
               zg=zpi-zp0(j)-zmap1(j1)
               xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
               yc=           co(2,2)*yg+co(2,3)*zg
               zc=                      co(3,3)*zg
               rsq=xc*xc+yc*yc+zc*zc
               auxa=1.d0/dsqrt(rsq)
               r6=rsq*rsq*rsq/sigma**6.
               gp=(apnd*alch_fact(jj,1) + r6)
               gp1=1.d0/gp
               gp2=gp1*gp1 
               gp0=(apnd*alch_fact0(jj,1) + r6)
               gp10=1.d0/gp0
               gp20=gp10*gp10 
               gr=6.d0*r6*auxa
               ssvir=alch_fact(jj,2)*epsmin*(2.d0*gp1*gp2 - gp2) 
               ssvirr=gr*ssvir
               auxc=epsmin*(gp2 - gp1)
               auxc0=epsmin*(gp20 - gp10)
               ucon=alch_fact(jj,2)*auxc
               ucon0=alch_fact0(jj,2)*auxc0
               cmap2(j1)=cmap2(j1)+ucon
               uconf(typeij)=uconf(typeij)+ucon*swrs(j1)
               uconf0(typeij)=uconf0(typeij)+ucon0*swrs(j1)
               qforce=ssvirr*swrs(j1)*auxa
c               if(i1.lt.9) write(6,101),i1,j,epsmin,sigma,dsqrt(rsq)
c     &              ,auxc*1.D4,epsmin*(2.d0*gp1*gp2 - gp2)*1.d4,
c     &              alch_fact(jj,2),alchemylj(jj),qforce,ssvirr,ssvir,gr
c101            format(2i5,2g15.5,f9.4, 3x,2G15.5,2f8.1,2g15.5)  
               auxc=chrgei*charge(j)*auxa
               ssvir=alch_factq(jj,2)*auxc
               ssvir0=alch_factq0(jj,2)*auxc
               cmap2(j1)=cmap2(j1)+ssvir
               ucoul(typeij)=ucoul(typeij)+ssvir*swrs(j1)
               ucoul0(typeij)=ucoul0(typeij)+ssvir0*swrs(j1)
               qforce=qforce+ssvir*swrs(j1)*auxa*auxa
               emvir=qforce
               fpx(i1)=fpx(i1)+qforce*xc
               fpy(i1)=fpy(i1)+qforce*yc
               fpz(i1)=fpz(i1)+qforce*zc
               fpx(j)=fpx(j)-qforce*xc
               fpy(j)=fpy(j)-qforce*yc
               fpz(j)=fpz(j)-qforce*zc

               st1 = st1+emvir*xc*xxcm
               st2 = st2+emvir*xc*yycm
               st3 = st3+emvir*xc*zzcm
               st4 = st4+emvir*yc*xxcm
               st5 = st5+emvir*yc*yycm
               st6 = st6+emvir*yc*zzcm
               st7 = st7+emvir*zc*xxcm
               st8 = st8+emvir*zc*yycm
               st9 = st9+emvir*zc*zzcm

            END DO
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
2002     noff=nnlpp(1+n)+1
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

      dwrk=0.d0 
c--   compute the work in the alchemical transofrmation
      do i=1,3
        auxa=ucoul(i)-ucoul0(i)
        auxb=uconf(i)-uconf0(i)
        dwrk=dwrk+auxa+auxb
      end do

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

      RETURN
      END
