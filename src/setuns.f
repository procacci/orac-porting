      SUBROUTINE setuns(time)

************************************************************************
*                                                                      *
*     Set up the unit system used in the program. Values are           *
*     assigned to some fundamental physical constants. All             *
*     these variables are left untouched throughout the program.       *
*                                                                      *
*--------------------- ARGUMENT ---------------------------------------*
*                                                                      *
*     TIME    :  Timestep size.                                        *
*                                                                      *
*----- Last update  04/08/89/ -----------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8   time

*==================== EXECUTABLE STATEMENTS ============================

      pi=4.0D0*DATAN(1.0D0)
      twopi=2.0D0*pi
      avogad=6.022140857d23
      boltz=1.38064852d-23
      gascon=8.3144598d0
      planck=6.62607004d-34
      elechg=1.60217662d-19
      epso=8.854187817d-12
      boxl=2.0D0
      unitm=1.0D0/(avogad*1000.0D0)
      unitl=1.0D-10
      unitt=1.D-15
      unite=unitm*(unitl/unitt)**2
      unitc=4.0D0*pi*epso*unitl*unite/(elechg*elechg)
      unitp=(unite/unitl**3)/1.0D6
      efact=unite*avogad
      ttstep=time
*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
