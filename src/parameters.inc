
      INTEGER untp,noato,corrp,backp,bendp,torsp,itorp,linkg
      INTEGER atpep,hydp,norsre,noatre,bondp
      PARAMETER(untp=200,noato=350,corrp=1,backp=10,hydp=20,bendp=10000,
     x          torsp=3000,itorp=350,atpep=m6,linkg=100,norsre=3,
     x          noatre=4,bondp=5000)

!=======================================================================
!
!     UNTP   :  Maximum number of residue units.
!     NOATO   :  Maximum number of atoms per residue.
!     CORRP   :  Maximum number of solvent-macromolecule
!                GR's per residue.
!     BACKP   :  Maximum number of backbone atoms per residue.
!     HYDP    :  Maximum number of hydrogen bond for residue.
!     BENDP   :  Maximum number of bending parameters.
!     BONDP   :  Maximum number of bond parameters.
!     TORSP   :  Maximum number of proper torsion parameters.
!     ITORP   :  Maximum number of improper torsion parameters.
!     ATPEP   :  Maximum number of types in the model potential.
!     LINKG   :  Maximum number of extra links.
!     NORSRE  :  Maximum number of rigid body per residue
!     NOATRE  :  Maximum number of secondary atoms per rigid bodies of
!                a residue
!
!=======================================================================

      CHARACTER*7 alpha(noato,untp),natype(noato,untp),jtor(4,noato,untp
     &     ),jitor(4,noato,untp),jbnd(2,noato,untp),cnat(3,untp)
     &     ,alphb(backp,untp),jacc(2,hydp,untp),jdon(2,hydp,untp)
     &     ,albnd(3,bendp),altor(4,torsp),alito(4,itorp),alhyd(2,atpep
     &     *(atpep+1)/2),xbond(2,linkg),xtor(4,linkg*3),xitor(4,linkg)
     &     ,jbend(3,noato,untp),stpr(4,norsre,untp),stsc(noatre,norsre
     &     ,untp),albon(2,bondp)

      INTEGER nctpe(noato,untp),natop(untp),nmbo(untp),ntor(untp)
     &     ,nitor(untp),nback(untp),nacc(untp),ndon(untp),ntor2(torsp)
     &     ,jngrp(untp),rbond(2,linkg),rtor(4,3*linkg),ritor(4,linkg)
     &     ,jgrppt(2,noato,untp),nstpr(norsre,untp),nstsc(norsre,untp)
     &     ,nbend(untp),nstrg(untp),xnbond,xntor,xnitor
      INTEGER lpbnd,lptor,lpito,lpbon
      REAL*8  qchge(noato,untp),qmss(noato,untp),pbnd1(bendp)
     &     ,massu(atpep),pbnd2(bendp),pbnd3(bendp),pbnd4(bendp)
     &     ,pbnd5(bendp),pbnd6(bendp),pbnd7(bendp),ptor1(torsp)
     &     ,pito1(itorp),pito2(itorp),pito3(itorp),pbnd(bendp,7)
     &     ,ptor(torsp,2),pito(itorp,3),potj(torsp,4),pbon1(bondp)
     &     ,pbon2(bondp),pbon(bondp,2)

!------------------- VARIABLES IN COMMON -------------------------------


      COMMON /dymam/ qchge,qmss,pbon1,pbon2,pbnd1,pbnd2,pbnd3,pbnd4
     &     ,pbnd5,pbnd6,pbnd7,ptor1,pito1,pito2,pito3,pbon,pbnd,potj
     &     ,ptor,pito,massu,lpbnd,lptor,lpito,lpbon,xnbond,xntor,xnitor
     &     ,nctpe,natop,nmbo,nbend,ntor,nitor,nback,nacc,ndon,ntor2
     &     ,rbond,rtor,ritor,jngrp,jgrppt,nstpr,nstsc,nstrg,alpha,natype
     &     ,jtor,jitor,jbnd,jbend,cnat,alphb,jacc,jdon,albnd,altor,alito
     &     ,albon,alhyd,stpr,stsc,xbond,xtor,xitor



