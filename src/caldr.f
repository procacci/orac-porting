      SUBROUTINE caldr(xyz0,xyzfit,wc,work,nato)
      IMPLICIT none
      REAL*8  xyz0(3,*),xyzfit(3,*),work(*),wc(*)
      LOGICAL near0
      INTEGER nato

      INTEGER n,m,na
      REAL*8  xd,yd,zd,dr,sum
      DATA m/0/

      sum=0.0D0
      DO n=1,nato
          IF(.NOT. near0(wc(n))) THEN
              xd=xyz0(1,n)-xyzfit(1,n)
              yd=xyz0(2,n)-xyzfit(2,n)
              zd=xyz0(3,n)-xyzfit(3,n)
              dr=xd*xd+yd*yd+zd*zd
              work(n)=dr
          ELSE
              work(n)=0.0D0
          END IF
      END DO

      RETURN
      END
