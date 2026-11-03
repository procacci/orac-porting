      SUBROUTINE chkewl(oc,rsqcut,alphal,rkcut,v)

************************************************************************
*                                                                      *
*                                                                      *
*     Check convergence of the Ewald summation and construct           *
*     vectors to be used on the sum on the resiprocal lattice.         *
*                                                                      *
*     OC      :  Transform the crordinates to simulation box      (I)  *
*                frame.                                                *
*                >> real*8 OC(3,3) <<                                  *
*     QCHG    :  List of atomic charges on the solvent molecules. (I)  *
*                >> real*8  QCHG(NTS) <<                               *
*     CHRGPW  :  List of atomic charges on the solute molecules   (I)  *
*                as seen by the solvent, packed by molecule.           *
*                >> real*8  CHRGPW(NTAP) <<                            *
*     NTS     :  Total number of sites on the solvent molecules.  (I)  *
*     NTAP    :  Total number of sites on the molecules of the    (I)  *
*                solute.                                          (I)  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     RKCUT   :  Ewald sum cutoff over the k-vectors.             (I)  *
*     V       :  Volume of the simulation box.                    (I)  *
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
*---- Last update 06/07/94 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CEA, Saclay France 1993                *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use unit
      use fourier

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  oc(3,3),alphal,rkcut,v,rsqcut

*-------------------- LOCAL VARIABLES ----------------------------------

      LOGICAL near0
      INTEGER i,l,m,n,mp,np,loff,moff,noff,ncount,
     x        ilmax,immax,inmax
      REAL*8  ro(3,3),rkscut,tworcl,tworcm,tworcn,rolm1,rolm2,rolm3,gl,
     x        gm,gn,qpefac,qfrfac,qa,qb,dunc,rl,rm,rn,rksq,
     x        kmax,emax,akm,eps,klmax,kmmax,knmax,a1,a2,a3,a4,a5,rspcut,
     x        qt,qp,alphar,expcst,erfcst,ucoul

      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*================ EXECUTABLE STATEMENTS ================================


      WRITE(kprint,10000)
!=======================================================================
!----- Check convergence on the direct lattice -------------------------
!=======================================================================

      rspcut=DSQRT(rsqcut)
      alphar=alphal*rspcut
      qt=1.0d0/(1.0d0+qp*alphar)
      expcst=DEXP(-alphar*alphar)
      erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
      ucoul=erfcst/rspcut
      WRITE(kprint,20000)
      WRITE(kprint,20100) alphal,rspcut,ucoul

!=======================================================================
!----- Check convergence on the reciprocal lattice ---------------------
!=======================================================================

      ncount=0
      klmax=0.0D0
      kmmax=0.0D0
      knmax=0.0D0
      ilmax=0
      immax=0
      inmax=0
      emax=1.0D10
      kmax=0.0D0
      dunc=DSQRT(unitc)
      qpefac=2.0d0/unitc
      qfrfac=-2.0d0
      qa=2.0D0*twopi/v
      qb=1.0d0/(4.0d0*alphal*alphal)
      rkscut=rkcut**2
      eps=rkscut/100.0D0
      tworcl=2.0d0*pi/boxl
      tworcm=2.0d0*pi/boxl
      tworcn=2.0d0*pi/boxl


      mmint=mofset
      nmint=nofset+1
      DO 10 i=1,3
          ro(i,1)=oc(1,i)
          ro(i,2)=oc(2,i)
          ro(i,3)=oc(3,i)
10    CONTINUE
      DO 60 loff=1,lmaxt
          l=loff-lofset
          rl=tworcl*DFLOAT(l)
          IF(loff.NE.1) THEN
              mmint=1
              nmint=1
          END IF
          DO 70 moff=mmint,mmaxt
              m=moff-mofset
              mp=IABS(m)+1
              rm=tworcm*dfloat(m)
              rolm1=ro(1,1)*rl+ro(1,2)*rm
              rolm2=ro(2,1)*rl+ro(2,2)*rm
              rolm3=ro(3,1)*rl+ro(3,2)*rm
              DO 90 noff=nmint,nmaxt
                  n=noff-nofset
                  np=iabs(n)+1
                  rn=tworcn*dfloat(n)
                  gl=rolm1+ro(1,3)*rn
                  gm=rolm2+ro(2,3)*rn
                  gn=rolm3+ro(3,3)*rn
                  rksq=gl*gl+gm*gm+gn*gn
                  IF(rksq .LE. rkscut) THEN
                      ncount=ncount+1
                      akm=qa*DEXP(-qb*rksq)/rksq
                      IF(emax .GT. akm) emax=akm
                      IF(kmax .LT. DSQRT(rksq)) kmax=DSQRT(rksq)
                      IF(knmax .LT. DABS(gn) ) THEN
                         knmax=DABS(gn)
                         inmax=n
                      END IF
                      IF(kmmax .LT. DABS(gm)) THEN
                         kmmax=DABS(gm)
                         immax=m
                      END IF
                      IF(klmax .LT. DABS(gl)) THEN
                         klmax=DABS(gl)
                         ilmax=l
                      END IF
                  END IF
90            CONTINUE
              nmint=1
70        CONTINUE
60    CONTINUE
      WRITE(kprint,30000)
      WRITE(kprint,30100) rkcut,LMAX, MMAX, NMAX,ilmax,immax,inmax,
     x                    klmax,kmmax,knmax,ncount,
     x                    emax
      WRITE(kprint,40000)

10000 FORMAT(//80('=')/'=',78x,'='
     x   /'=                        Convergence',
     x' Test on Ewald Sum                         ='/
     x       '=',78x,'='/80('=')//)
20000 FORMAT(/16x,'**********************************'/
     x        16x,'*                                *'/
     x        16x,'*           Direct Sum           *'/
     x        16x,'*                                *'/
     x        16x,'**********************************'//)

20100 FORMAT(1x,'  Alpha             = ',f10.3,'      1/Ang '/
     x       1x,'  R_c               = ',f10.3,'      Ang '/
     x       1x,'  Energy at Cutoff  = ',e12.5,
     x          '      Ang.-1    '//)

30000 FORMAT(/16x,'**********************************'/
     x        16x,'*                                *'/
     x        16x,'*         Reciprocal Sum         *'/
     x        16x,'*                                *'/
     x        16x,'**********************************'//)

30100 FORMAT(1x,'  Rkcut                            =',f10.2,
     x      20x,'   1/Ang '/
     x       1x,'  Set k-vectors LMAX, MMAX, NMAX   =',3(4x,i6)/
     x       1x,'  Smallest Possible Choice         =',3(4x,i6)/
     x       1x,'  Largest k-vectors along l, m, n  =',3f10.5,
     x          '   1/Ang '/
     x       1x,'  Total Number of K-vectors        =',4x,i6/
     x       1x,'  Energy at Cutoff                 =',1x,e12.4,
     x      20x,'   Ang-1'/)

40000 FORMAT(//80('=')/80('=')//)
      

*================ END EXECUTABLE STATEMENTS ============================

      RETURN
      END
