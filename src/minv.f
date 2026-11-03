      SUBROUTINE minv(n,a,d,l,m)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION a(*),l(*),m(*)
c
c        ...............................................................
c
c        if a real version of this routine is desired, the
c        c in column 1 should be removed from the real
c        statement which follows.
c
c     REAL a,d,biga,hold
c
c        the c must also be removed from real statements
c        appearing in other routines used in conjunction with this
c        routine.
c
      IF(n.EQ.1) go to 1999
c
c        search for largest element
c
      d=1.0
      nk=-n
      DO 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      DO 20 j=k,n
      iz=n*(j-1)
      DO 20 i=k,n
      ij=iz+i
   10 IF(abs(biga)-abs(a(ij))) 15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 CONTINUE
c
c        interchange rows
c
      j=l(k)
      IF(j-k) 35,35,25
   25 ki=k-n
      DO 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=m(k)
      IF(i-k) 45,45,38
   38 jp=n*(i-1)
      DO 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by minus pivot (value of pivot element is
c        contained in biga)
c
   45 IF(biga) 48,46,48
   46 d=0.0
      RETURN
   48 DO 55 i=1,n
      IF(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 CONTINUE
c
c        reduce matrix
c
      DO 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      DO 65 j=1,n
      ij=ij+n
      IF(i-k) 60,65,60
   60 IF(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 CONTINUE
c
c        divide row by pivot
c
      kj=k-n
      DO 75 j=1,n
      kj=kj+n
      IF(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 CONTINUE
c
c        product of pivots
c
      d=d*biga
c
c        replace pivot by reciprocal
c
      a(kk)=1.0/biga
   80 CONTINUE
c
c        final row and column interchange
c
      k=n
  100 k=(k-1)
      IF(k) 150,150,105
  105 i=l(k)
      IF(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      DO 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      IF(j-k) 100,100,125
  125 ki=k-n
      DO 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 RETURN
 1999 a(1)=1./a(1)
      RETURN
      END
