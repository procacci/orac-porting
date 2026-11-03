      FUNCTION brent(ax,bx,cx,f,tol,xmin,xp0,yp0,zp0
     &     ,fpx,fpy,fpz)

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      EXTERNAL f
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
      REAL*8  ax,bx,cx,tol,xmin,f,brent

*----------------------- LOCAL VARIABLES -------------------------------

      INTEGER itmax,iter
      REAL*8  fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      REAL*8  a,b,d,e,etemp,cgold,zeps
      parameter (itmax=100,cgold=.3819660,zeps=1.0e-10)
  

*======================= EXECUTABLE STATEMENTS =========================


      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.0D0
      fx=f(x,xp0,yp0,zp0,fpx,fpy,fpz)
      fv=fx
      fw=fx
      DO iter=1,itmax
          xm=0.5*(a+b)
          tol1=tol*DABS(x)+zeps
          tol2=2.*tol1
          IF(DABS(x-xm).le.(tol2-.5*(b-a))) goto 3
          IF(DABS(e).gt.tol1) then
              r=(x-w)*(fx-fv)
              q=(x-v)*(fx-fw)
              p=(x-v)*q-(x-w)*r
              q=2.*(q-r)
              IF(q.gt.0.) p=-p
              q=DABS(q)
              etemp=e
              e=d
              IF(DABS(p).ge.DABS(.5*q*etemp).or.p.le.q*(a-x).or. 
     *            p.ge.q*(b-x)) goto 1
              d=p/q
              u=x+d
              IF(u-a.lt.tol2 .or. b-u.lt.tol2) d=DSIGN(tol1,xm-x)
              goto 2
          END IF
1         IF(x.ge.xm) then
              e=a-x
          ELSE
              e=b-x
          END IF
          d=cgold*e
2         IF(DABS(d).ge.tol1) then
              u=x+d
          ELSE
              u=x+DSIGN(tol1,d)
          END IF
          fu=f(u,xp0,yp0,zp0,fpx,fpy,fpz)
          IF(fu.le.fx) then
              IF(u.ge.x) then
                  a=x
              ELSE
                  b=x
              END IF
              v=w
              fv=fw
              w=x
              fw=fx
              x=u
              fx=fu
          ELSE
              IF(u.lt.x) then
                  a=u
              ELSE
                  b=u
              END IF
              IF(fu.le.fw .or. w.eq.x) then
                  v=w
                  fv=fw
                  w=u
                  fw=fu
              ELSE IF(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                  v=u
                  fv=fu
              END IF
          END IF
      END DO
      write(6,*)"Errormsg from NUMREC: brent exceed maximum iterations."
      STOP 
3     xmin=x 
      brent=fx 

      RETURN
      END
