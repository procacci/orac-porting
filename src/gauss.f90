      SUBROUTINE gauss(tmp,mass,nato,nmol,vx,vy,vz)

!***********************************************************************
!                                                                      *
!                                                                      *
!     GAUSS will generate a set of gaussianly distributed              *
!     random velocities with average zero and standard deviation       *
!     DSQRT(tmp).                                                      *
!                                                                      *
!                                                                      *
!     TMP     :  Set temperature.                                 (I)  *
!     OC      :  Transform the coordinates to simulation          (I)  *
!                box frame.                                            *
!                >> real*8 OC(3,3) <<                                  *
!     MASS    :  List of atomic masses.                           (I)  *
!                >> real*8 MASS(NATO*NMOL) <<                          *
!     NATO    :  Number of atoms in the molecule with non         (I)  *
!                zero mass.                                            *
!     NCNST   :  Number of constraints per molecule.              (I)  *
!     NMOL    :  Number of molecules.                             (I)  *
!     VX      :  List of atomic velocities.                       (O)  *
!     VY         >> real*8 VX(NATO*NMOL), ... <<                       *
!     VZ                                                               *
!                                                                      *
! --- Last update 02/11/13 --------------------------------------------*
!                                                                      *
!     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
!                                                                      *
!     EXTERNAL RANF.                                                   *
!                                                                      *
!************************************************************************

!======================= DECLARATIONS ==================================

      use unit
#ifdef _MPI_
      use rem
      use orac_mpi
#endif

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,nmol
      REAL*8  mass(*),vx(*),vy(*),vz(*),tmp

!-------------------- LOCAL VARIABLES ----------------------------------

      EXTERNAL ranf,duni
      INTEGER i,nf,nts
      REAL*8  u1,u2,u3,u4,sig,dummy,t1,t2,t3,ranf,dustar,tvel

!-------------------- VARIABLES IN COMMONS -----------------------------

      LOGICAL near0
      DATA dummy/0.0D+0/

!==================== EXECUTABLE STATEMENTS ============================

      nts=nato*nmol
#ifdef _MPI_
!     seed is different for each trajectory and identical for each iproc_f 
      dummy=dustar(iproc)  
!     write(kprint,*) "Iproc from gauss", iproc,dummy,ranf() 
#else
      dummy=dustar(1)
#endif
      DO 10 i=1,nts
          tvel=0.0D0
          IF(.NOT. near0(mass(i))) tvel=boltz*tmp/(unite*mass(i))
          sig=DSQRT(tvel)
          u1=ranf(dummy)
          u2=ranf(dummy)
          u3=ranf(dummy)
          u4=ranf(dummy)
          IF(u1 .EQ. 0.0D0) THEN
              u1=1.0D-15
          END IF
          IF(u2 .EQ. 0.0D0) THEN
              u2=1.0D-15
          END IF
          IF(u3 .EQ. 0.0D0) THEN
              u3=1.0D-15
          END IF
          IF(u4 .EQ. 0.0D0) THEN
              u4=1.0D-15
          END IF
          t1=DSQRT(-2.0d0*DLOG(u1))*DCOS(2.0d0*pi*u2)
          t2=DSQRT(-2.0d0*DLOG(u1))*DSIN(2.0d0*pi*u2)
          t3=DSQRT(-2.0d0*DLOG(u3))*DCOS(2.0d0*pi*u4)
          vx(i)=t1*sig
          vy(i)=t2*sig
          vz(i)=t3*sig
10    CONTINUE

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the kinetic energy,
! which can be used to rescale the velocities.
! The procedure can be applied to all atoms or to smaller groups.
! If it is applied to intersecting groups in sequence, the kinetic energy
! that is given as an input (kk) has to be up-to-date with respect to the previous rescalings.
!
! When applied to the entire system, and when performing standard molecular dynamics (fixed c.o.m. (center of mass))
! the degrees of freedom of the c.o.m. have to be discarded in the calculation of ndeg,
! and the c.o.m. momentum HAS TO BE SET TO ZERO.
! When applied to subgroups, one can chose to:
! (a) calculate the subgroup kinetic energy in the usual reference frame, and count the c.o.m. in ndeg
! (b) calculate the subgroup kinetic energy with respect to its c.o.m. motion, discard the c.o.m. in ndeg
!     and apply the rescale factor with respect to the subgroup c.o.m. velocity.
! They should be almost equivalent.
! If the subgroups are expected to move one respect to the other, the choice (b) should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces to an istantaneous
! randomization of the kinetic energy, as described in paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator against discretization errors.
! The easiest recipe is:
! htilde = h + conint
! where h is the total energy (kinetic + potential)
! and conint is a quantity accumulated along the trajectory as minus the sum of all the increments of kinetic
! energy due to the thermostat.
!
function resamplekin(kk,sigma,ndeg,taut)
  implicit none
  real*8               :: resamplekin
  real*8,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  real*8,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  real*8,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  real*8 :: factor,rr
  real*8, external :: gasdev
  factor=exp(-1.0/taut)
  rr = gasdev()
  resamplekin = kk + (1.0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk) &
               + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor)
!  write(103,123) taut,factor,kk,sigma,rr,ndeg,resamplekin,sumnoises(ndeg-1),sigma*(sumnoises(ndeg-1)+rr**2)/ndeg,(sumnoises(ndeg-1)+rr**2)/ndeg
!123 format(2g15.3,2f15.5,f10.5,i10,4f15.5)
contains 

double precision function sumnoises(nn)
  implicit none
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
  real*8, external :: gamdev,gasdev
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0*gamdev(nn/2)
  else
    sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

end function resamplekin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING ROUTINES ARE TRANSCRIBED FROM NUMERICAL RECIPES

double precision function gamdev(ia)
! gamma-distributed random number, implemented as described in numerical recipes

implicit none
integer, intent(in) :: ia
integer j
real*8 am,e,s,v1,v2,x,y
real*8, external :: ran1
if(ia.lt.1)pause 'bad argument in gamdev'
if(ia.lt.6)then
  x=1.
  do 11 j=1,ia
    x=x*ran1()
11  continue
  x=-log(x)
else
1 v1=2.*ran1()-1.
    v2=2.*ran1()-1.
  if(v1**2+v2**2.gt.1.)goto 1
    y=v2/v1
    am=ia-1
    s=sqrt(2.*am+1.)
    x=s*y+am
  if(x.le.0.)goto 1
    e=(1.+y**2)*exp(am*log(x/am)-s*y)
  if(ran1().gt.e)goto 1
endif
gamdev=x
end

double precision function gasdev()
! gaussian-distributed random number, implemented as described in numerical recipes

implicit none
integer, save :: iset = 0
real*8, save :: gset
real*8, external :: ran1
real*8 fac,rsq,v1,v2
if(iset==0) then
1       v1=2.*ran1()-1.0d0
  v2=2.*ran1()-1.0d0
  rsq=v1**2+v2**2
  if(rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
end if
end

FUNCTION ran1()
! random number generator
INTEGER IA,IM,IQ,IR,NTAB,NDIV
REAL*8 ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
INTEGER, SAVE :: idum=0 !! ATTENTION: THE SEED IS HARDCODED
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
  do 11 j=NTAB+8,1,-1
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11  continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
      END

      
      subroutine removep(ntap,vpx,vpy,vpz,xmp,massinfty)
        
        !---     remove net linear momentum  
        implicit NONE

        real*8 vpx(*),vpy(*),vpz(*),xmp(*),massinfty  !args
        integer ntap                                  !arg
        integer l                                     !local
        real*8 vbx,vby,vbz,xmtot,bgg                      !local  
        
        vbx = 0.d0
        vby = 0.d0
        vbz = 0.d0
        xmtot=0.d0
        do l = 1, ntap
           if(xmp(l).gt.massinfty) THEN 
              vpx(l)=0.d0
              vpy(l)=0.d0
              vpz(l)=0.d0
           END IF
           vbx = vbx + vpx(l)*xmp(l)
           vby = vby + vpy(l)*xmp(l)
           vbz = vbz + vpz(l)*xmp(l)
           xmtot=xmp(l)+xmtot
        end do
        bgg=1.d0/(xmtot)
        vbx = vbx*bgg
        vby = vby*bgg
        vbz = vbz*bgg
        
        do l=1,ntap
           if(xmp(l).lt.massinfty) THEN
              vpx(l) =  vpx(l)-vbx
              vpy(l) =  vpy(l)-vby
              vpz(l) =  vpz(l)-vbz
           end if
        end do
        
    end subroutine removep

!=================================================================================
    subroutine bussi_thermostat(vpx,vpy,vpz,vcax,vcay,vcaz,temp, & 
         efact,gascon,mass,massp,massinfty,ebussi,friction,ntap,nprot,co,cnstpp) 
!=================================================================================
!    Bussi thermostat. To be called at the end of the outermost time step
!    (e.g. right before mtsmd_dump.inc [call to mts_test routine] mtsmd_dump
!    calls kinetic(Vs) that passes ucek and pucek to mts_times. ebussi is the 
!    only variable mts_test must know to implement HTilde and this must be passed
!    to mts_test from mtsmd_dump.
!    N.B: ATOMIC VELOCITIES vpx,vpy,vpz, are in the ORTHOGNAL system with zero COM velocity 
!         for proteins/molecules
!         COM VELOCITIES vcax,vcay,vcaz are in the CRYSTAL system. 
!=================================================================================
!---------------------Declarations------------------------------------------------
      implicit NONE 
      
!ARGS
      REAL*8  :: vpx(*),vpy(*),vpz(*),mass(*),massp(*),vcax(*),vcay(*),vcaz(*)
      REAL*8  :: massinfty,ebussi,friction,co(3,3),gascon,temp,efact
      INTEGER :: ntap,nprot,cnstpp
!LOCAL 
     REAL*8 :: vxyz(3),ek,tto,ttcm,ekcm,scale,sigma,ekt,resamplekin,rtemp
     INTEGER :: ndeg,i,cnstinfty

     REAL*8 vcmx,vcmxa,masst
!---------------------executable statements--------------------------------

     ekt=0.5*gascon*temp  ! 1/2 K_B T (j mol-1)  (target temperature)
     tto=0.d0


!--> ! compute current values of mol-frame kinetic energy
     cnstinfty=0
     DO i=1,ntap
        if(mass(i).GT.massinfty) cnstinfty=cnstinfty+3
        tto=tto+0.50*mass(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
     end do
     tto=tto*efact
!    scale relative (0 CM) velocities within proteins/molecules
     ndeg = 3*ntap-3*nprot - cnstpp -cnstinfty
     sigma=ndeg*ekt  
     ek=resamplekin(tto,sigma,ndeg,friction)  
     scale=dsqrt(ek/tto)
     rtemp=2.d0*ek/ndeg/gascon
     do i=1,ntap   
        vpx(i)=vpx(i)*scale
        vpy(i)=vpy(i)*scale
        vpz(i)=vpz(i)*scale   ! uniform scale factor preserve 0 CM velocities and constraints 
     end do
     ebussi=ebussi-(ek - tto )  ! bussi term (in Kj mol-1) for HTilde
!     write(105,*) "bussi intra", ebussi,-(ek - tto ) 
!     write(101,123) ebussi,scale,ek,tto,rtemp,friction
!123  FORMAT(5F15.3,G15.5) 

!--->! compute current values of COM kinetic energy of proteins/molecule
     ttcm=0.d0
     DO i=1,nprot
        vxyz(1)=co(1,1)*vcax(i)+co(1,2)*vcay(i)+co(1,3)*vcaz(i)
        vxyz(2)=co(2,1)*vcax(i)+co(2,2)*vcay(i)+co(2,3)*vcaz(i)
        vxyz(3)=co(3,1)*vcax(i)+co(3,2)*vcay(i)+co(3,3)*vcaz(i)
        ttcm=ttcm+0.5D0*massp(i)*(vxyz(1)**2+vxyz(2)**2+vxyz(3)**2)
     END DO
     ttcm=ttcm*efact
!    scale COM velocities of the proteins/molecules
     ndeg = 3*nprot-3 
     sigma=ndeg*ekt
     ekcm=resamplekin(ttcm,sigma,ndeg,friction)  
     scale=dsqrt(ekcm/ttcm)
     vcmx=0.d0
     vcmxa=0.d0
     masst=0.d0
     call removep(nprot,vcax,vcay,vcaz,massp,massinfty)
     do i=1,nprot
        vcax(i)=vcax(i)*scale
        vcay(i)=vcay(i)*scale
        vcaz(i)=vcaz(i)*scale   ! uniform scale factor can be safely
                                ! applied to crystal frame velocities
        vcmx=vcmx+massp(i)*vcax(i)
        masst=masst+massp(i)
        vcmxa=vcmxa +abs(vcax(i))
     end do
!     write(102,64) vcmx/masst,vcmxa/nprot,masst
!64   FORMAT(5G15.5)
     ebussi=ebussi-(ekcm - ttcm )  ! bussi term for HTilde
     rtemp=2.d0*ekcm/ndeg/gascon
!     write(102,123) ebussi,scale,ekcm,ttcm,rtemp,friction
   end subroutine bussi_thermostat
