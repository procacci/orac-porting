      SUBROUTINE anneal(fact,vx,vy,vz,nato)

      IMPLICIT none

      REAL*8 fact,vx(*),vy(*),vz(*)
      INTEGER  nato

      INTEGER i

      DO i=1,nato
         vx(i)=vx(i)*fact
         vy(i)=vy(i)*fact
         vz(i)=vz(i)*fact
      END DO
      
      RETURN
      END

