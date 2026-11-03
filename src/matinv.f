      SUBROUTINE matinv(n,nfirst,a,ainv,det)
      REAL*8 a(nfirst,n),ainv(nfirst,n)
      REAL*8 det
      INTEGER lt(21),mt(21)
      LOGICAL near0
      DO 30 i=1,n
      DO 30 j=1,n
      ainv(j,i)=a(j,i)
 30   CONTINUE
      CALL minv(n,ainv,det,lt,mt)
      IF(.not.near0(det))RETURN
      print*,' warning the matrix is singular '
      det=0.0d0
      DO 20 i=1,n
      DO 20 j=1,n
      ainv(j,i)=0.0d0
 20   CONTINUE
      RETURN
      END
