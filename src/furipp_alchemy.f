CCCCC  GFS CCCC        #include "pressure.H"
      SUBROUTINE furipp_alchemy(ss_index,oc,xp0,yp0,zp0,charge,charge0
     &     ,ntap,atomp,grppt,alphal,rkcut,v,urecp_slt,urecp_slv,urecp_ss
     &     ,xpcm,ypcm,zpcm,fppx,fppy,fppz,co,stresstc,dwrk)

************************************************************************
*                                                                      *
*                                                                      *
*     The contributions from the sum over the reciprocal lattice       *
*     to the solute-solute electrostatic energy and forces             *
*     are computed by FURIPP. The energy is evaluated by               *
*     taking the square of the geometric structure factor              *
*     multiplied by the sites charges. Addition and subtraction        *
*     formulas for cosines and sines are used.                         *
*                                                                      *
*                                                                      *
*     OC      :  Transform the crordinates to simulation box      (I)  *
*                frame.                                                *
*                >> real*8 OC(3,3) <<                                  *
*     XP0     :  Coordinates for the solvent molecules, packed    (I)  *
*     YP0        by molecule.                                          *
*     ZP0        >> real*8  XP0(NTAP), ... <<                          *
*     CHARGE  :  List of atomic charges on the solute molecules   (I)  *
*                as seen by the solvent, packed by molecule.           *
*                >> real*8  CHARGE(NTAP) <<                            *
*     NTAP    :  Total number of sites on the molecules of the    (I)  *
*                solute.                                          (I)  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     RKCUT   :  Ewald sum cutoff over the k-vectors.             (I)  *
*     V       :  Volume of the simulation box.                    (I)  *
*     URECP   :  Contribution from the reciprocal lattice sum     (O)  *
*                to the electrostatic energy.                          *
*     VIRMM   :  Contribution from the reciprocal lattice sum     (O)  *
*                to the virial.                                        *
*     STRESSTC(3,3) : Coulomb contribution to the stress tensor        *
*     EMVIR         :  Virial                                          *
*                                                                      *
*                                                                      *
*=================== Definition of the sum boundaries =================*
*                                                                      *
*     The sum over the reciprocal space is carried out in the 3        *
*     dimensions from -LMAX, 0, LMAX  -MMAX, 0, MMAX  -NMAX, 0, NMAX   *
*                                                                      *
*     However only half of the reciprocal lattice vectors are used     *
*     The constant lofset, mofset, nofset, lmaxt, mmaxt and nmaxt      *
*     are defined from LMAX, MMAX and NMAX.                            *
*                                                                      *
*     lofset=1                                                         *
*     mofset = MMAX + 1                                                *
*     nofset = NMAX + 1                                                *
*     lmaxt  = LMAX + 1                                                *
*     mmaxt  = 2*MMAX + 1                                              *
*     nmaxt  = 2*NMAX + 1                                              *
*                                                                      *
*========================== WORK ARRAYS ===============================*
*                                                                      *
*     CLM     :  >> real*8  CLM(NA) <<                                 *
*     SLM     :  >> real*8  SLM(NA) <<                                 *
*     EXPIKR  :  >> real*8  EXPIKR(NA) <<                              *
*     EXPIKI  :  >> real*8  EXPIKI(NA) <<                              *
*     ELR     :  >> real*8  ELR(NA,LMAX+1) <<                          *
*     ELI     :  >> real*8  ELI(NA,LMAX+1) <<                          *
*     EMR     :  >> real*8  EMR(NA,MMAX+1) <<                          *
*     EMI     :  >> real*8  EMI(NA,MMAX+1) <<                          *
*     ENR     :  >> real*8  ENR(NA,NMAX+1) <<                          *
*     ENI     :  >> real*8  ENI(NA,NMAX+1) <<                          *
*     NA .GE. NTAP                                                     *
*     If this is not satisfied the routine aborts the run.             *
*                                                                      *
*---- Last update 01/07/93 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1993               *
*                                                                      *
*     EXTERNALS XERROR.                                                *
*                                                                      *
************************************************************************

C======================= DECLARATIONS ==================================

      use unit
      use parst
      use fourier

      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,ss_index(*)
      REAL*8  oc(3,3),urecp_slt,urecp_slv,urecp_ss,alphal,rkcut,v
      REAL*8  xp0(*),yp0(*),zp0(*),charge(*),charge0(*),fppx(*),fppy(*)
     &     ,fppz(*),dwrk

C     Quantities needed to compute the stress tensor
      REAL*8  stresstc(3,3),co(3,3)
      REAL*8  xpcm(*),ypcm(*),zpcm(*)


C-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nd
      REAL*8 zero,one,two,three,six,twelve
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,
     x          three=3.0D0,six=6.0D0,twelve=12.0D0)
      PARAMETER (nd=m1)
      REAL*8  clm(nd),slm(nd),expikr(nd),expiki(nd),elr(nd,LMAX+1)
     &     ,eli(nd,LMAX+1),emr(nd,MMAX+1),emi(nd,MMAX+1),enr(nd,NMAX+1)
     &     ,eni(nd,NMAX+1),fpx(m1),fpy(m1),fpz(m1)
      COMMON / rag1 / clm,slm,expikr,expiki,elr,eli,emr,emi,enr,eni,fpx
     &     ,fpy,fpz

C-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER i,j,l,m,n,mp,np,loff,moff,noff,lmn,lmn1,typei
      REAL*8  ro(3,3),qpefac,qfrfac,qa,qb,rkscut,tworcl,tworcm,tworcn,qx
     &     ,qy,qz,rl,rm,rn,rolm1,rolm2,rolm3,gl,gm,gn,dsnm,emii,enii
     &     ,rksq,fl,fn,fm,rk,akm,dsnn,enet(3),qc,qforce,sumrp(2),sumip(2
     &     ),urecp(3),sumipt,sumrpt,enett,ene0,sumr0,sumi0

      REAL*8  sgx(m1),sgy(m1),sgz(m1),sx,sy,sz,aux1,aux2,aux3
      INTEGER k,grppt(2,*),atomp(*),count,ia,mia
      REAL*8 stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9






C================ EXECUTABLE STATEMENTS ================================

C=======================================================================
C---- Check the dimensions of the work arrays --------------------------
C=======================================================================

      IF(nd.LT.ntap) THEN
          errmsg=' IN FURIPP : Dimensions of the work arrays are insuffi
     xcient. ABORT CC '
          WRITE(6,*) errmsg
          STOP
      END IF

      do i=1,ntap
        fpx(i)=0.0d0
        fpy(i)=0.0d0
        fpz(i)=0.0d0
      end do

      IF(DABS(rkcut) .LT. 1.0D-2) RETURN
      qpefac=two
      qfrfac=-two
      qa=qpefac*twopi/v
      qb=one/(4.0D0*alphal*alphal)
      rkscut=rkcut**2
      tworcl=two*pi/boxl
      tworcm=two*pi/boxl
      tworcn=two*pi/boxl
      mmint=mofset
      nmint=nofset+1
#ifndef PRESSURE
      DO i=1,3
         urecp(i)=0.0D+0
      END DO

#else /* PRESSURE */
      ene0=0.d0
      urecp=zero
      stc1=zero
      stc2=zero
      stc3=zero
      stc4=zero
      stc5=zero
      stc6=zero
      stc7=zero
      stc8=zero
      stc9=zero

C     Computes atomic positions with respect to the centre of mass in 
C     cartesian coordinates 

      DO i=1,ntap
         j=atomp(i)
         sx=xp0(i)-xpcm(j)
         sy=yp0(i)-ypcm(j)
         sz=zp0(i)-zpcm(j)
         sgsx(i)=sx*co(1,1)+sy*co(1,2)+sz*co(1,3)
         sgsy(i)=sx*co(2,1)+sy*co(2,2)+sz*co(2,3)
         sgsz(i)=sx*co(3,1)+sy*co(3,2)+sz*co(3,3)
      END DO

#endif
      DO 10 i=1,3
          ro(i,1)=oc(1,i)
          ro(i,2)=oc(2,i)
          ro(i,3)=oc(3,i)
10    CONTINUE
      DO 20 i=1,ntap
          qx=tworcl*xp0(i)
          qy=tworcm*yp0(i)
          qz=tworcn*zp0(i)
          elr(i,1)=one
          emr(i,1)=one
          enr(i,1)=one
          eli(i,1)=zero
          emi(i,1)=zero
          eni(i,1)=zero
          elr(i,2)=DCOS(qx)
          emr(i,2)=DCOS(qy)
          enr(i,2)=DCOS(qz)
          eli(i,2)=DSIN(qx)
          emi(i,2)=DSIN(qy)
          eni(i,2)=DSIN(qz)
20    CONTINUE
      DO lmn=2,lmax
          lmn1=lmn+1
          DO i=1,ntap
              elr(i,lmn1)=elr(i,lmn)*elr(i,2)-eli(i,lmn)*eli(i,2)
              eli(i,lmn1)=elr(i,lmn)*eli(i,2)+eli(i,lmn)*elr(i,2)
          END DO
      END DO
      DO lmn=2,mmax
          lmn1=lmn+1
          DO i=1,ntap
              emr(i,lmn1)=emr(i,lmn)*emr(i,2)-emi(i,lmn)*emi(i,2)
              emi(i,lmn1)=emr(i,lmn)*emi(i,2)+emi(i,lmn)*emr(i,2)
          END DO
      END DO
      DO lmn=2,nmax
          lmn1=lmn+1
          DO i=1,ntap
              enr(i,lmn1)=enr(i,lmn)*enr(i,2)-eni(i,lmn)*eni(i,2)
              eni(i,lmn1)=enr(i,lmn)*eni(i,2)+eni(i,lmn)*enr(i,2)
          END DO
      END DO
      DO 50 loff=1,lmaxt
          l=loff-lofset
          rl=tworcl*dfloat(l)
          IF(loff.NE.1) THEN
              mmint=1
              nmint=1
          END IF
          DO 60 moff=mmint,mmaxt
              m=moff-mofset
              mp=IABS(m)+1
              rm=tworcm*dfloat(m)
              rolm1=ro(1,1)*rl+ro(1,2)*rm
              rolm2=ro(2,1)*rl+ro(2,2)*rm
              rolm3=ro(3,1)*rl+ro(3,2)*rm
              dsnm=DBLE(SIGN(1,m))
              DO 70 i=1,ntap
                  emii=dsnm*emi(i,mp)
                  clm(i)=elr(i,loff)*emr(i,mp)-eli(i,loff)*emii
                  slm(i)=elr(i,loff)*emii     +eli(i,loff)*emr(i,mp)
70            CONTINUE
              DO 80 noff=nmint,nmaxt
                  n=noff-nofset
                  np=IABS(n)+1
                  rn=tworcn*dfloat(n)
                  gl=rolm1+ro(1,3)*rn
                  gm=rolm2+ro(2,3)*rn
                  gn=rolm3+ro(3,3)*rn
                  rksq = gl*gl + gm*gm + gn*gn        
                  IF(rksq.LE.rkscut) THEN
                      fl=oc(1,1)*gl+oc(1,2)*gm+oc(1,3)*gn
                      fm=oc(2,1)*gl+oc(2,2)*gm+oc(2,3)*gn
                      fn=oc(3,1)*gl+oc(3,2)*gm+oc(3,3)*gn
                      rk=DSQRT(rksq)
                      akm=qa*DEXP(-qb*rksq)/rksq
#ifndef PRESSURE
                      DO i=1,2
                         sumrp(i)=zero
                         sumip(i)=zero
                      END DO
                      dsnn=dble(sign(1,n))
#else /* PRESSURE */
                      aux1=two*(one+qb*rksq)/rksq
                      DO i=1,2
                         sumrp(i)=zero
                         sumip(i)=zero
                      END DO
                      dsnn=DBLE(SIGN(1,n))
#endif 
                      
                      sumr0=0.d0
                      sumi0=0.d0
                      DO 90 i=1,ntap
                         typei=ss_index(i)
                         enii=dsnn*eni(i,np)
                         expikr(i)=enr(i,np)*clm(i)-enii*slm(i)
                         expiki(i)=enii*clm(i)    +enr(i,np)*slm(i)
                         sumrp(typei)=sumrp(typei)+charge(i)*expikr(i)
                         sumip(typei)=sumip(typei)+charge(i)*expiki(i)
                         sumr0=sumr0+charge0(i)*expikr(i)
                         sumi0=sumi0+charge0(i)*expiki(i)
90                    CONTINUE
                      enet(1)=akm*(sumrp(1)*sumrp(1)+sumip(1)*sumip(1))
                      enet(2)=2.0D0*akm*(sumrp(2)*sumrp(1)+sumip(2)
     &                     *sumip(1))
                      enet(3)=akm*(sumrp(2)*sumrp(2)+sumip(2)*sumip(2))
                      ene0=akm*(sumr0*sumr0 + sumi0*sumi0)
                      sumipt=sumip(1)+sumip(2)
                      sumrpt=sumrp(1)+sumrp(2)
                      enett=0.D0
                      DO i=1,3
                         urecp(i)=urecp(i)+enet(i)
                         enett=enett+enet(i)
                      END DO
#ifndef PRESSURE
                      DO 100 j=1,ntap
                          qc=qfrfac*akm*charge(j)
                          qforce=qc*(expikr(j)*sumipt-expiki(j)*sumrpt)
                          fpx(j)=fpx(j)+fl*qforce
                          fpy(j)=fpy(j)+fm*qforce
                          fpz(j)=fpz(j)+fn*qforce
#else /* PRESSURE */
c      Coulombic contribution to the stress tensor
                      stc1=stc1+enett*(one-aux1*gl*gl) 
                      stc2=stc2+enett*( -aux1*gl*gm) 
                      stc3=stc3+enett*( -aux1*gl*gn) 
                      stc4=stc4+enett*( -aux1*gm*gl) 
                      stc5=stc5+enett*(one-aux1*gm*gm) 
                      stc6=stc6+enett*( -aux1*gm*gn) 
                      stc7=stc7+enett*( -aux1*gn*gl) 
                      stc8=stc8+enett*( -aux1*gn*gm) 
                      stc9=stc9+enett*(one-aux1*gn*gn) 

                      aux1=zero
                      aux2=zero
                      aux3=zero

                      DO 100 k=1,ntap
                          qc=qfrfac*akm*charge(k)
                          qforce=qc*(expikr(k)*sumipt-expiki(k)*sumrpt)
                          fpx(k)=fpx(k)+fl*qforce
                          fpy(k)=fpy(k)+fm*qforce
                          fpz(k)=fpz(k)+fn*qforce
                          aux1=aux1+qforce*sgx(k)
                          aux2=aux2+qforce*sgy(k)
                          aux3=aux3+qforce*sgz(k)
#endif 
100                   CONTINUE
#ifdef PRESSURE
                      stc1=stc1-aux1*gl
                      stc2=stc2-aux2*gl
                      stc3=stc3-aux3*gl
                      stc4=stc4-aux1*gm
                      stc5=stc5-aux2*gm
                      stc6=stc6-aux3*gm
                      stc7=stc7-aux1*gn
                      stc8=stc8-aux2*gn
                      stc9=stc9-aux3*gn

#endif 
                  END IF
80            CONTINUE
              nmint=1
60        CONTINUE
50    CONTINUE
#ifdef PRESSURE

      stresstc(1,1)=stresstc(1,1)+stc1
      stresstc(1,2)=stresstc(1,2)+stc2
      stresstc(1,3)=stresstc(1,3)+stc3
      stresstc(2,1)=stresstc(2,1)+stc4
      stresstc(2,2)=stresstc(2,2)+stc5
      stresstc(2,3)=stresstc(2,3)+stc6
      stresstc(3,1)=stresstc(3,1)+stc7
      stresstc(3,2)=stresstc(3,2)+stc8
      stresstc(3,3)=stresstc(3,3)+stc9
#endif 

      urecp_slt=urecp(1)
      urecp_ss=urecp(2)
      urecp_slv=urecp(3)
      dwrk=urecp(1)+urecp(2)+urecp(2)-ene0

      DO i=1,ntap
         fppx(i)=fppx(i)+co(1,1)*fpx(i)+co(1,2)*fpy(i)+co(1,3)*fpz(i)
         fppy(i)=fppy(i)+co(2,1)*fpx(i)+co(2,2)*fpy(i)+co(2,3)*fpz(i)
         fppz(i)=fppz(i)+co(3,1)*fpx(i)+co(3,2)*fpy(i)+co(3,3)*fpz(i)
      END DO

C================ END EXECUTABLE STATEMENTS ============================

      RETURN
      END
