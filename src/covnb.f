      SUBROUTINE covnb(pot1,pot2,pot3,pot4,pota,potb,npot,iz,hyd1,hyd2,
     x     lhyd,mass,tmass,ntap,pot6,pot12,pot146,pot1412,c6,c12
     &     ,type_table,lj_fudge,n1)

************************************************************************
*                                                                      *
*     The non-bonded interaction parameters for the macromolecule      *
*     are converted to program units. Mixing rules are used for        *
*     non hydrogen bond interactions. Additionally, the                *
*     macromolecule total mass is computed.                            *
*                                                                      *
*                                                                      *
*     POT1    :  List of alpha parameters.                    (INPUT)  *
*                >> real*8 POT1(*) <<                                  *
*     POT2    :  List of N-i parameters.                      (INPUT)  *
*                >> real*8 POT2(*) <<                                  *
*     POT3    :  List of R-i parameters.                      (INPUT)  *
*                >> real*8 POT3(*) <<                                  *
*     NPOT    :  Number of potential parameter in the         (INPUT)  *
*                model potential.                                      *
*     IZ      :  Mixing rules flag. IZ=0 no mixing rules      (INPUT)  *
*                are used. IZ=1 mixing rules are used.                 *
*     HYD1    :  List of hydrogen bond parameters.    (INPUT/OUTPUT)   *
*                >> real*8 HYD1(*) <<                                  *
*     HYD2    :  List of hydrogen bond parameters.    (INPUT/OUTPUT)   *
*                >> real*8 HYD2(*) <<                                  *
*     LHYD    :  Number of hydrogen bond parameters.                   *
*     MASS    :  List of atomic masses.               (INPUT/OUTPUT)   *
*                >> real*8   MASS(NTAP) <<                             *
*     TMASS   :  Total mass of the macromolecule.            (OUTPUT)  *
*     NTAP    :  Number of atoms making up the                (INPUT)  *
*                macromolecule.                                        *
*     POT12   :  List of repulsive L-J parameters in         (OUTPUT)  *
*                program units.                                        *
*                >> real*8 POT12(N1) <<                                *
*     POT6    :  List of repulsive L-J parameters in         (OUTPUT)  *
*                program units.                                        *
*                >> real*8 POT6(N1) <<                                 *
*     N1      :  Physical dimension of POT6/POT12             (INPUT)  *
*                Larger or equal to NPOT*(NPOT+1)/2                    *
*                                                                      *
*---- Last update 04/11/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n1
      INTEGER lhyd,ntap,npot,iz,type_table(n1,*)
      REAL*8 pot1(*),pot2(*),pot3(*),pot4(*),pot12(*),pot6(*),pota(*)
     &     ,potb(*),c6(*),c12(*),pot146(*),pot1412(*)
      REAL*8 hyd1(*),hyd2(*),mass(*),tmass,lj_fudge,aux

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,ij,naux
      REAL*8 tmpi,tmpj,emin,rmin,qener,fact,tmpj1,tmpj2
      LOGICAL near0
      CHARACTER*80 errmsg
      EXTERNAL near0

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*--------- Change the potential parameters units to program units ------
*--------- Apply mixing rules ------------------------------------------
*=======================================================================

      fact=1.0D0/(2.0D0**(1.0D0/6.0D0))
      aux =1.d-12
      naux=0
      IF(iz.EQ.0) THEN
         DO 10 i=1,npot
            DO 20 j=i,npot
               ij=j*(j-1)/2 + i

*=======================================================================
*----   Regular L-J parameters
*=======================================================================

               qener=DSQRT(DABS(pot2(i)*pot2(j)))
               tmpi=4.0D0*1000.0D0*qener*4.184/(unite*avogad)
               if(ffield_type(1:5).eq."amber") THEN 
                 tmpj=pot1(i)+pot1(j)
                 tmpj=(tmpj*fact)**6
                 pot6(ij)=tmpi*tmpj
                 pot12(ij)=tmpi*tmpj**2
               else if(ffield_type(1:4).eq."opls") THEN 
                 tmpj1=(2*pot1(i)*fact)**3
                 tmpj2=(2*pot1(j)*fact)**3
                 pot6(ij)=tmpi*tmpj1*tmpj2
                 pot12(ij)=tmpi*(tmpj1*tmpj2)**2
               end if
               type_table(i,j)=ij
               type_table(j,i)=ij

*=======================================================================
*----   Jorgensen parameters type
*=======================================================================

               qener=DSQRT(DABS(potb(i)*potb(j)))
               tmpi=4.0D0*1000.0D0*qener*4.184/(unite*avogad)
               tmpj=pota(i)+pota(j)
               tmpj=(tmpj*fact)**6
               c6(ij)=tmpi*tmpj
               c12(ij)=tmpi*tmpj**2

*=======================================================================
*----   1-4 Lennard-Jones parameters
*=======================================================================

               qener=DSQRT(DABS(pot4(i)*pot4(j)))
               tmpi=4.0D0*1000.0D0*qener*4.184/(unite*avogad)
               if(ffield_type(1:5).eq."amber") THEN 
                 tmpj=pot3(i)+pot3(j)
                 tmpj=(tmpj*fact)**6
                 pot146(ij)=tmpi*tmpj
                 pot1412(ij)=tmpi*tmpj**2
               else if(ffield_type(1:5).eq."opls") THEN 
                 tmpj1=(2*pot3(i)*fact)**3
                 tmpj2=(2*pot3(j)*fact)**3
                 pot146(ij)=tmpi*tmpj1*tmpj2
                 pot1412(ij)=tmpi*(tmpj1*tmpj2)**2
               end if
               IF(abs(pot1412(ij)) .LT. aux) THEN 
                 pot146(ij)=pot6(ij)*lj_fudge
                 pot1412(ij)=pot12(ij)*lj_fudge
               END IF
20          CONTINUE
10       CONTINUE
      ELSE IF(iz.EQ.1) THEN

*=======================================================================
*--------- Change the potential parameters units to program units ------
*--------- Assume that mixing rules are not to be used. Convert only ---
*--------- POT12 and POT6 ----------------------------------------------
*=======================================================================

         DO i=1,npot
            DO j=i,npot
               ij=j*(j-1)/2 + i
               type_table(i,j)=ij
               type_table(j,i)=ij
               pot6(ij)=1000.0D0*pot6(ij)*4.184/(unite*avogad)
               pot12(ij)=1000.0D0*pot12(ij)*4.184/(unite*avogad)
*
*             MODIFIED BY MATTEO 23/06/96
*
               pot1412(ij)=pot12(ij)*lj_fudge
               pot146(ij)=pot6(ij)*lj_fudge
*
            END DO
         END DO
      END IF

*=======================================================================
*-------- Hydrogen bond parameters -------------------------------------
*=======================================================================

      IF(hyd2(1).LT.5.0D0) THEN
          DO 60 i=1,lhyd
              emin=hyd1(i)
              rmin=hyd2(i)
              hyd1(i)=-2.0D0*emin*rmin**6
              hyd2(i)=-3.0D0*emin*rmin**4
              hyd1(i)=1000.0D0*hyd1(i)*4.184/(unite*avogad)
              hyd2(i)=1000.0D0*hyd2(i)*4.184/(unite*avogad)
60        CONTINUE
      ELSE
          DO 70 i=1,lhyd
              hyd1(i)=1000.0D0*hyd1(i)*4.184/(unite*avogad)
              hyd2(i)=1000.0D0*hyd2(i)*4.184/(unite*avogad)
70        CONTINUE
      END IF

*=======================================================================
*-------- compute total mass -------------------------------------------
*=======================================================================

      tmass=0.0D0
      DO 80 i=1,ntap
          tmass=tmass+mass(i)
80    CONTINUE
      
      lj_fudge=1.0D0

      
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
