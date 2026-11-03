*===============================================================================
      SUBROUTINE scale_masses(ntap,mass,atomtoscale,velscale,nscalemass
     &     ,vpx,vpy,vpz)
*===============================================================================
c--   Scale atomic masses when starting form restart if required 
*===============================================================================


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  vpx(*),vpy(*),vpz(*),velscale(*),mass(*)
      integer  atomtoscale(*)
      INTEGER nscalemass,ntap

*----------------------- LOCAL VARIABLES ------------------------------

      INTEGER i,j

*==================== EXECUTABLE STATEMENTS ============================


      DO j=1,nscalemass
        i=atomtoscale(j)
        vpx(i)=vpx(i)*velscale(j)
        vpy(i)=vpy(i)*velscale(j)
        vpz(i)=vpz(i)*velscale(j)
        mass(i)=mass(i)/(velscale(j)**2)
      END DO

      RETURN 
      END
