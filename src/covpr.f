      SUBROUTINE covpr(volume,t,taut,taup,pext,mpr,wpr,compressibility)
      use unit
      IMPLICIT none

      REAL*8  mpr,wpr,t,volume,compressibility
      REAL*8  taut,taup,pext

      EXTERNAL near0
      LOGICAL near0

      REAL*8  c,omega,bulkmod,v


      taut=taut/unitt
      taup=taup/unitt
      pext=pext/unitp

      c=2.997925D10

C     Converts wavelenght in omega=1/T

      omega=twopi*c*wpr

C     Converts volume in SI

      v=volume*unitl**3

C     This is the bulk modulus of a perfect gas. 
C     We should find a better approximation for b-m.
C     From the expression bulkmod=pressure valid for  perfect gas
C     we obtain a BM = 1E6 dyne/cm^2 = 0.1 MPa = 1 atm 
C                 BM = 1.7E10 dyne/cm^2 for liquid CF4
C                 BM = 2.7E10 dyne/cm^2 for solid Ar
C                 BM = 2.4E11 dyne/cm^2 for solid NaCl 
C                 BM = 2.0E10 dyne/cm^2 for spc 300K (berendsen)
C       1GPa = 1E10 dyne/cm^2
C     bulkmod should be given in Pa. 
C     From simulation we obtain:
C      Substance     T/K       B/Pa
C     |-----------|--------|--------|
C         Ar        100      2.3E8   ;   4.47E8 (10000st)
C         H2O(spc)  300      2.7E9   


      bulkmod=1.0D6/compressibility

      IF (.NOT.near0(omega)) THEN
          mpr=12.0D0*bulkmod*v**(1.0D0/3.0D0)/(omega**2)
          mpr=mpr/unitm
      END IF

      IF (near0(omega)) THEN
          mpr=1.0D199
      END IF

      RETURN
      END
