      FUNCTION chebev(a,b,c,m,x)
      IMPLICIT none
      INTEGER m
      REAL*8  a,b,x,c(m),chebev
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.0D0) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.0D0*x-a-b)/(b-a)
      y2=2.0D0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
