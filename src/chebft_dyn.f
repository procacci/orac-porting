      SUBROUTINE chebft_dyn(c,n,fc)
      IMPLICIT none
      INTEGER n,NMAX
      REAL*8 c(n),fc(*),PI
      PARAMETER (NMAX=50, PI=3.141592653589793d0)
      INTEGER j,k
      REAL*8 fac,y,sum

      fac=2.0D0/DFLOAT(n)
      do 13 j=1,n
         sum=0.d0
         do 12 k=1,n
            sum=sum+fc(k)*cos((PI*(j-1))*((k-0.5d0)/n))
12       continue
         c(j)=fac*sum
13    continue
      return
      END
