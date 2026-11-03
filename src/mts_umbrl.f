      SUBROUTINE mts_umbrl(unbias,flag,iflag,xp0,yp0,zp0,uconf,nbone
     &     ,mback,nato,enout,omega,alpha,gr,fpx,fpy,fpz)

************************************************************************
*                                                                      *
*     Umbrella subroutine.                                             *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*---- Last update 03/10/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UCB, Berkeley 1991                     *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use parst

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      
      REAL*8  uconf,enout,omega,alpha,xp0(*),yp0(*),zp0(*)
      REAL*8  fpx(*),fpy(*),fpz(*)
      REAL*8  gr
      INTEGER nato,nbone,iflag,mback(*)
      LOGICAL flag,unbias

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,n,m,nbtot,initial
      REAL*8  upper,lower,coef,grp,fnbtot
      REAL*8  xc,yc,zc,xr,yr,zr,xpi,ypi,zpi,rsq
      REAL*8  fppx,fppy,fppz
      LOGICAL near0
      COMMON /rag1/ fppx(m1),fppy(m1),fppz(m1)
      DATA initial/0/
      SAVE initial

*-------------------- LOCAL VARIABLES ----------------------------------

      uconf=0.0D0
      gr=0.0D0
      nbtot=nbone*(nbone-1)/2
      fnbtot=1.0D0/DFLOAT(nbtot)

      DO i=1,nato
          fppx(i)=0.0D0
          fppy(i)=0.0D0
          fppz(i)=0.0D0
      END DO
      DO m=1,nbone
          i=mback(m)
          xpi=xp0(i)
          ypi=yp0(i)
          zpi=zp0(i)
          DO n=m+1,nbone
              j=mback(n)
              xr=xpi-xp0(j)
              yr=ypi-yp0(j)
              zr=zpi-zp0(j)
              rsq=xr**2+yr**2+zr**2
              gr=gr+rsq
              fppx(i)=fppx(i)-xr
              fppy(i)=fppy(i)-yr
              fppz(i)=fppz(i)-zr
              fppx(j)=fppx(j)+xr
              fppy(j)=fppy(j)+yr
              fppz(j)=fppz(j)+zr
          END DO
      END DO
      gr=gr*fnbtot
      grp=DSQRT(gr)
      IF(near0(enout)) THEN
         enout=grp
      END IF
      IF(initial .EQ. 0) THEN
         WRITE(*,1000) enout
      END IF

      IF(flag .OR. unbias) THEN
         IF(.NOT. unbias) THEN
            IF(iflag .LT. 0) THEN
            
*=======================================================================
*--- If the solute must be folded --------------------------------------
*=======================================================================

               IF(enout .GE. grp) THEN
                  enout=grp
               ELSE
                  coef=2.0D0*alpha*(enout-grp)*fnbtot/grp
                  uconf=alpha*(enout-grp)**2
                  DO m=1,nbone
                     i=mback(m)
                     fpx(i)=fpx(i)-coef*fppx(i)
                     fpy(i)=fpy(i)-coef*fppy(i)
                     fpz(i)=fpz(i)-coef*fppz(i)
                  END DO
               END IF

*=======================================================================
*--- If the solute must be unfolded ------------------------------------
*=======================================================================

            ELSE
               IF(enout .LE. grp) THEN
                  enout=grp
               ELSE
                  coef=2.0D0*alpha*(enout-grp)*fnbtot/grp
                  uconf=alpha*(enout-grp)**2
                  DO m=1,nbone
                     i=mback(m)
                     fpx(i)=fpx(i)-coef*fppx(i)
                     fpy(i)=fpy(i)-coef*fppy(i)
                     fpz(i)=fpz(i)-coef*fppz(i)
                  END DO
               END IF
            END IF
         ELSE
            coef=2.0D0*alpha*(enout-grp)*fnbtot/grp
            uconf=alpha*(enout-grp)**2
            DO m=1,nbone
               i=mback(m)
               fpx(i)=fpx(i)-coef*fppx(i)
               fpy(i)=fpy(i)-coef*fppy(i)
               fpz(i)=fpz(i)-coef*fppz(i)
            END DO
         END IF
      END IF
      gr=grp
      initial=1
       
!================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(//'----------- R_g = ', f12.4,' -------------'//)
      RETURN
      END
