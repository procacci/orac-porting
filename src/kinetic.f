      SUBROUTINE kinetic(ss_index,co,nato_slv,nmol,cnstpp_slv,ntap
     &     ,cnstpp,vpx,vpy,vpz,vpxcm,vpycm,vpzcm,nprot,protl,mass,massp
     &     ,wmtp,ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk
     &     ,cpress,isostress,vco,masspp,ucepr,temppra,massinfty)

************************************************************************
*                                                                      *
*   KINETIC computes all Kinetic energies and temperatures of both     *
*   solvent and solute, c.o.m velocities and rms of solvent molecules. *
*   KINETIC is called only when NRESPA > 1                             *
*                                                                      *
*     SOLVEN  :  Solvent flag.                                    (I)  *
*     PROTEI  :  Macromolecule flag.                              (I)  *
*     OC      :  Transform the coordinates to simulation box      (I)  *
*                frame.                                                *
*                >> real*8  OC(3,3)  <<                                *
*     CO      :  Transform the coordinates to orthogonal          (I)  *
*                crystallographic frame.                               *
*                >> real*8  OC(3,3)  <<                                *
*     XP0     :  Solute coordinates, packed by molecule.          (I)  *
*     YP0        >> real*8  XP0(*), ... <<                             *
*     ZP0                                                              *
*     VPX     :  Solute velocities , packed by molecule.          (I)  *
*     VPY        >> real*8  XP0(*), ... <<                             *
*     VPZ                                                              *
*     VCX     :  Solvent velocities, packed by sites.             (I)  *
*     VCY        >> real*8  XP0(*), ... <<                             *
*     VCZ                                                              *
*     WMSS    :  List of the solvent atomic masses, packed by     (I)  *
*                site.                                                 *
*                >> real*8  MASSA(*) <<                                *
*     WTMOL   :  Total mass of the solvent molecule.              (I)  *
*     MASS    :  List of solute atomic masses, packed by molecule.(I)  *
*                >> real*8  MASSB(*) <<                                *
*     WMTP    :  Total mass of the solute molecules.              (I)  *
*     NATOW   :  Number of atoms per solvent molecule.            (I)  *
*     NMOL    :  Number of solvent molecules.                     (I)  *
*     NTAP    :  Number of atomic sites on the solute.            (I)  *
*                                                                      *
*     XSMOVE  :  Solvent rms                                     (I/O) *
*     YSMOVE  :  >> real*8  XSMOVE(*), ... <<                    (I/O) *
*     ZSMOVE  :                                                  (I/O) *
*                                                                      *
*     UCEK    :  Kinetic energy of the solvent molecules.         (O)  *
*     PUCEK   :  Kinetic energy of the solute.                    (O)  *
*     TEMP    :  Total temperature of the system                  (O)  *
*     TEMPT   :  Translational temperature of the solvent.        (O)  *
*     TEMPR   :  Rotational + intra temperature of the solvent.   (O)  *
*     TEMPPR  :  Temperature of the solute                        (O)  *
*     TCM     :  Translational temperature of the solute.         (O)  *
*     RCM     :  Rotational temperature of the solute.            (O)  *
*     V0X     :  Solute c.o.m. velocities                         (O)  *
*     V0Y        >> real*8  XP0(*), ... <<                        (0)  *
*     V0Z                                                         (0)  *
*                                                                      *
*     NA .GE. NMOL*NATOW, NB .GE. NTAP, NC .GE. NMOL                   *
*     If this conditions are not satisfied the program aborts.         *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst
      use omp_integr, only:omp_timing,time_in,tcpu_in,typei

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER   nato_slv,nmol,ntap,nprot,protl(*),ss_index(*),cnstpp_slv
     &     ,cnstpp
      REAL*8    mass(*),wmtp,co(3,3),vco(3,3),ucepr,masspp(*),temppra
     &     ,massinfty
      REAL*8    vpx(*),vpy(*),vpz(*)
      REAL*8    ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk(3,3)
      REAL*8  vpxcm(*),vpycm(*),vpzcm(*),massp(*)
      LOGICAL cpress,isostress

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,nf,nfslv,nfslt,type,iii,jjj,count,j
      REAL*8  vxyz(3)
      REAL*8  tto(2),ttcm(2),tucek,mtot
      INTEGER na,nb,offset,cnstinfty(2)
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      PARAMETER (na=npm,nb=m1)

*==================== EXECUTABLE STATEMENTS ============================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
      
      DO i=1,2
         tto(i)=0.d0
         cnstinfty(i)=0
      END DO
#ifdef _OMP_
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(type) 
!$OMP&  REDUCTION(+:tto,cnstinfty)
#endif
      DO i=1,ntap
         type=ss_index(i)
         tto(type)=tto(type)+(0.50D+0)*mass(i)*(vpx(i)**2+vpy(i)**2
     &        +vpz(i)**2)
         if(mass(i).gt.massinfty) THEN 
           cnstinfty(type)=cnstinfty(type)+3
         end if
      end do   
#ifdef _OMP_
!$OMP END PARALLEL DO
#endif

      DO iii=1,3
         DO jjj=1,3
            stresstk(iii,jjj) = 0.0D0
         END DO
      END DO
      DO i=1,2
         ttcm(i)=0.0D0
      END DO
      count=0
#ifdef _OMP_
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(type,vxyz) 
!$OMP& REDUCTION(+:ttcm,stresstk)
      DO i=1,nprot
         type=typei(i)
         vxyz(1)=co(1,1)*vpxcm(i)+co(1,2)*vpycm(i)+co(1,3)*vpzcm(i)
         vxyz(2)=co(2,1)*vpxcm(i)+co(2,2)*vpycm(i)+co(2,3)*vpzcm(i)
         vxyz(3)=co(3,1)*vpxcm(i)+co(3,2)*vpycm(i)+co(3,3)*vpzcm(i)
         ttcm(type)=ttcm(type)+0.5D0*massp(i)*(vxyz(1)**2+vxyz(2)**2
     &        +vxyz(3)**2)
         DO iii=1,3
            DO jjj=1,3
              if(massp(i).lt.massinfty) 
     &             stresstk(iii,jjj) = stresstk(iii,jjj) + 
     x                            massp(i)*vxyz(iii)*vxyz(jjj)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
#else
      DO i=1,nprot
         offset=protl(count+1)
         type=ss_index(protl(count+1+1))
         vxyz(1)=co(1,1)*vpxcm(i)+co(1,2)*vpycm(i)+co(1,3)*vpzcm(i)
         vxyz(2)=co(2,1)*vpxcm(i)+co(2,2)*vpycm(i)+co(2,3)*vpzcm(i)
         vxyz(3)=co(3,1)*vpxcm(i)+co(3,2)*vpycm(i)+co(3,3)*vpzcm(i)
         ttcm(type)=ttcm(type)+0.5D0*massp(i)*(vxyz(1)**2+vxyz(2)**2
     &        +vxyz(3)**2)
         DO iii=1,3
            DO jjj=1,3
              if(massp(i).lt.massinfty) 
     &             stresstk(iii,jjj) = stresstk(iii,jjj) + 
     x                            massp(i)*vxyz(iii)*vxyz(jjj)
            END DO
         END DO
         count=count+offset+1
      END DO
#endif

      nf=3*ntap - cnstpp - cnstinfty(1) - cnstinfty(2)
      tucek=(tto(1)+ttcm(1)+tto(2)+ttcm(2))*efact
      pucek=(tto(1)+ttcm(1))*efact  
      ucek=(tto(2)+ttcm(2))*efact
      nfslv=3*nato_slv*nmol-cnstpp_slv - cnstinfty(2)
      nfslt=nf-nfslv


      temp=2.0D0*tucek/(gascon*nf)
      IF(nfslt .NE. 0) THEN
         temppr=2.0D0*pucek/(gascon*nfslt)
      ELSE
         temppr=0.0D0
      END IF

      tcm=ttcm(1)+ttcm(2)
      tcm=2.0D0*tcm*efact/(gascon*3.0D0*DFLOAT(nprot))
      IF(nf-3*nprot .NE. 0) THEN
         rcm=2.0D0*(tto(1)+tto(2))*efact/(gascon*DFLOAT(nf-3*nprot))
      ELSE
         rcm=0.0D0
      END IF
      tempt=tcm
      tempr=rcm

      IF(cpress) THEN
         ucepr=0.0D0
         DO i=1,3
            DO j=i,3
               ucepr=ucepr+masspp(i)*vco(i,j)*vco(i,j)
            END DO
         END DO
         ucepr=0.5D0*ucepr*efact
         IF(isostress) THEN
            temppra=2.0D0*ucepr/gascon
         ELSE
            temppra=(2.0D0/6.0D0)*ucepr/gascon
         END IF
      END IF

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF
      RETURN
      END
