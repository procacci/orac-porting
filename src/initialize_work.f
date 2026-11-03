!----------------------------------------------------------------------------
      subroutine initialize_work(xp0,yp0,zp0,xpcm,ypcm,zpcm,fpx
     &       ,fpy,fpz,wrk0) 
!----------------------------------------------------------------------------
!- Computes starting work when Ewald is on. The reference state is the
!  isolated molecule in the starting PBC box (if constant pressure =
!  true). So it must include all intrasolute electrostatics form BOTH
!  direct and reciprocal space. Such reference state has a (small)
!  inter-images electfrostatic interacion due to even multipoles.  The
!  starting state in an alchemical trasformation for a molecule to be
!  switched and the final state for a molecule to be switched off is the
!  following hybrid state: full intrasolute charges (direct and
!  reciprocal), inter-images solute-solute electrostatics but also a
!  spurious solute-solvent recirpocal term that must be subtraced out.
!  such term arises as one cannot set the charges of th solute to zero
!  for the reciprocal space without impacting on the intrasolute
!  electrostatics. This routine computes the starting/on or ending/off
!  work as follows:  
 
! The hydrib state has the correct solvent-solvent and solute-solute 
! electrostatic energies but an excess solute-solvent reciprocal contribution 
! the starting work for an "on" transformation is hence given by
! DWRK0=Energy_of_the_hybrid_state -True_energy of the reference system 
! where
! Energy_of_reference state = E(intrasolute)_dir + E(intrasolute)_recp +
! E(Intrasolute)_self  +  E(intrasolvent)_dir + E(intrasolvent)_recp +
! E(Intrasolvent__self) 
! Energy_of_hybrid state = E(intrasolute)_dir + E(intrasolute)_recp +
! E(intrasolute)_self + E(intrasolvent_dir) + E(intrasolvent_recp) +
! + E(intrasolvent)_selft + E(solute-solvent)recp  
!                           ^^^^^^^^^^^^^^^^^^^^^  
! This routine computes the true_energy (recp+self) of the reference state
! The energy of  the hybrid state (recp+self) is computed in mtsmd
!----------------------------------------------------------------------------
      use cpropar
      use parst
      use spme

      implicit none 
      real*8 wrk0,chrge0(m1),wrk0t,wrk0g,wrkref,wrk0s
      REAL*8  xp0(*),yp0(*),zp0(*),xpcm(*),ypcm(*),zpcm(*),fpx(*),fpy(*)
     &     ,fpz(*),self_slt,self_slv,eer,urcp,urcs,urcsp,virial(3,3)
      integer  i,nat
c---  compute reciprocal+self energy of the reference system

*===========================================================================================
c---- COMPUTE RECP+SELF OF GHOST/SOLUTE ATOMS (lambda=1): set the
c---- charges of the true atoms to zero. N.B. subsystem to be switched
c---- on may be smaller than the ORAC solute
*==========================================================================================
      do i=1,ntap 
        chrge0(i)= dabs(lambdaq(i))*chrge(i)
      end do

c--   recp
      eer=0.d0
      urcs=0.d0
      urcp=0.d0
      urcsp=0.d0
      IF(.NOT. pme) THEN
        CALL furipp(ss_index,oc,xp0,yp0,zp0,chrge0,ntap,atomp,grppt
     &       ,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm,zpcm,fpx
     &       ,fpy,fpz,co,virial)
        wrk0g=urcp+urcs+urcsp
      ELSE
        call fft_pme(ntap,xp0,yp0,zp0,xpcm,ypcm,zpcm,chrge0,co,oc
     &       ,volume,alphal,pme_order,nfft1,nfft2,nfft3,eer,fpx,fpy
     &       ,fpz,virial,atomp,grppt,pressure,rkcut)
        wrk0g=eer 
      END IF 


*=======================================================================
c---- COMPUTE RECP+SELF OF TRUE/SOLVENT ATOMS (lambda=0):set the charges of
c---- ghost atoms to zero.
*=======================================================================

      do i=1,ntap 
        chrge0(i)= (1-dabs(lambdaq(i)))*chrge(i)
      end do

c--   recp
      eer=0.d0
      urcs=0.d0
      urcp=0.d0
      urcsp=0.d0
      IF(.NOT. pme) THEN
         CALL furipp(ss_index,oc,xp0,yp0,zp0,chrge0,ntap,atomp,grppt
     &        ,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm,zpcm,fpx
     &        ,fpy,fpz,co,virial)
         wrk0t=urcp+urcs+urcsp
      ELSE
         call fft_pme(ntap,xp0,yp0,zp0,xpcm,ypcm,zpcm,chrge0,co,oc
     &        ,volume,alphal,pme_order,nfft1,nfft2,nfft3,eer,fpx,fpy
     &        ,fpz,virial,atomp,grppt,pressure,rkcut)
         wrk0t=eer 
      END IF 
      wrkref=wrk0t+wrk0g

*=======================================================================
c---- COMPUTE RECP+SELF OF TRUE+ghots SOLVENT+solvent ATOMS (all lambdas=0)
*=======================================================================

      do i=1,ntap 
        chrge0(i)= chrge(i)
      end do

c--   recp
      eer=0.d0
      urcs=0.d0
      urcp=0.d0
      urcsp=0.d0
      IF(.NOT. pme) THEN
         CALL furipp(ss_index,oc,xp0,yp0,zp0,chrge0,ntap,atomp,grppt
     &        ,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm,zpcm,fpx
     &        ,fpy,fpz,co,virial)
         wrk0s=urcp+urcs+urcsp
      ELSE
         call fft_pme(ntap,xp0,yp0,zp0,xpcm,ypcm,zpcm,chrge0,co,oc
     &        ,volume,alphal,pme_order,nfft1,nfft2,nfft3,eer,fpx,fpy
     &        ,fpz,virial,atomp,grppt,pressure,rkcut)
         wrk0s=eer 
      END IF 
!      write(kprint,100) wrk0t*1.d4,wrk0g*1.d4,wrk0s*1.d4
!100   format("solvent-rcp",f10.4," solute-rcp",f10.4, 
!     &     " solute-slv-recp",f10.4) 
      wrk0=wrk0s-wrkref  ! ref-> hybrid work (on process);   
      RETURN

      END 
