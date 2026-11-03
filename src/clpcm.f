      SUBROUTINE clpcm(mass,tmass,xa0,ya0,za0,nsite,ncm,xcm,ycm,zcm)

************************************************************************
*                                                                      *
*     Calculate the centre of mass for a given number of molecules.    *
*                                                                      *
*     MASS    :  List of atomic masses for the molecule.      (INPUT)  *
*                >> real*8 MASS(NSITE) <<                              *
*     TMASS   :  Total molecular mass.                        (INPUT)  *
*     XA0     :  Coordinates of the atoms for the set of      (INPUT)  *
*     YA0        molecules considered.                                 *
*     ZA0        >> real*8 XA0(NSITE*NCM), YA0(NSITE*NCM),             *
*                          ZA0(NSITE*NCM)  <<                          *
*     NSITE   :  Number of atomic sites on the molecule.      (INPUT)  *
*     NCM     :  Number of molecules.                         (INPUT)  *
*     XCM     :  Coordinates of the centres of mass.         (OUTPUT)  *
*     YCM        >> real*8 XCM(NCM), YCM(NCM), ZCM(NCM) <<             *
*     ZCM                                                              *
*                                                                      *
*----- Last update 04/20/89 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nsite,ncm
      REAL*8 tmass
      REAL*8 mass(*),xa0(*),ya0(*),za0(*),xcm(*),ycm(*),zcm(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,m,n
      REAL*8 fact

*==================== EXECUTABLE STATEMENTS ============================

      tmass=0.0D0
      DO i=1,nsite
         tmass=tmass+mass(i)
      END DO

      DO 10 m=1,ncm
          xcm(m)=0.0d0
          ycm(m)=0.0d0
          zcm(m)=0.0d0
10    CONTINUE
      n=0
      DO 20 i=1,nsite
          fact=mass(i)/tmass
          DO 30 m=1,ncm
              xcm(m)=xcm(m)+fact*xa0(n+m)
              ycm(m)=ycm(m)+fact*ya0(n+m)
              zcm(m)=zcm(m)+fact*za0(n+m)
30        CONTINUE
          n=n+ncm
20    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
