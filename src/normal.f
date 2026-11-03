      SUBROUTINE normal(a,n)
      IMPLICIT none
      REAL*8 a(*)
      LOGICAL near0
      INTEGER n
 
      REAL*8 aux
      INTEGER i

      aux=0.0D0
      DO i=1,n
          aux=aux+a(i)
      END DO
      IF(near0(aux)) aux=1.0D0
      DO i=1,n
          a(i)=a(i)/aux
      END DO
      RETURN
      END
