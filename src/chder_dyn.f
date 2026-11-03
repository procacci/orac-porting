      SUBROUTINE chder_dyn(a,b,c,cder,n)
      IMPLICIT none
      INTEGER n
      REAL*8  a,b,c(n),cder(n)
      INTEGER j
      REAL*8  con

      cder(n)=0.0D0
      cder(n-1)=2*(n-1)*c(n)
      do 11 j=n-2,1,-1
        cder(j)=cder(j+2)+2*j*c(j+1)
11    continue
      con=2.0D0/(b-a)
      do 12 j=1,n
        cder(j)=cder(j)*con
12    continue
      return
      END
