      SUBROUTINE comp_abmd_fdiss(iflag,enout,diss_list,alpha,xp0,yp0
     &     ,zp0,co,fpx,fpy,fpz,uconf,rsp)

************************************************************************
*   Time-stamp: <97/04/29 10:21:26 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iflag,diss_list(2,*)
      REAL*8  alpha,enout,xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
     &     ,co(3,3),rsp,uconf

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER k,i
      LOGICAL near0
      REAL*8  xd,yd,zd,xd1,yd1,zd1,xd2,yd2,zd2,xc,yc,zc,rsq,coef,aux1
     &     ,aux2
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*---- Compute distance -------------------------------------------------
*=======================================================================

      uconf=0.0D0
      aux1=1.0D0/DFLOAT(diss_list(1,1))
      aux2=1.0D0/DFLOAT(diss_list(2,1))
      xd1=0.0D0
      yd1=0.0D0
      zd1=0.0D0
      DO k=1,diss_list(1,1)
         i=diss_list(1,k+1)
         xd1=xd1+xp0(i)
         yd1=yd1+yp0(i)
         zd1=zd1+zp0(i)
      END DO
      xd1=xd1*aux1
      yd1=yd1*aux1
      zd1=zd1*aux1
      xd2=0.0D0
      yd2=0.0D0
      zd2=0.0D0
      DO k=1,diss_list(2,1)
         i=diss_list(2,k+1)
         xd2=xd2+xp0(i)
         yd2=yd2+yp0(i)
         zd2=zd2+zp0(i)
      END DO
      xd2=xd2*aux2
      yd2=yd2*aux2
      zd2=zd2*aux2

      xd=xd1-xd2
      yd=yd1-yd2
      zd=zd1-zd2
      xd=xd-2.0D0*PBC(xd)
      yd=yd-2.0D0*PBC(yd)
      zd=zd-2.0D0*PBC(zd)

      xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
      yc=co(2,1)*xd+co(2,2)*yd+co(2,3)*zd
      zc=co(3,1)*xd+co(3,2)*yd+co(3,3)*zd
      rsq=xc**2+yc**2+zc**2
      rsp=DSQRT(rsq)
      IF(near0(enout)) enout=rsp

      IF(iflag .LT. 0) THEN
            
*=======================================================================
*--- If the solute must be folded --------------------------------------
*=======================================================================

         IF(enout .GE. rsp) THEN
            enout=rsp
         ELSE
            coef=-2.0D0*alpha*(rsp-enout)/rsp
            uconf=alpha*(rsp-enout)**2
            DO k=1,diss_list(1,1)
               i=diss_list(1,k+1)
               fpx(i)=fpx(i)+coef*xc*aux1
               fpy(i)=fpy(i)+coef*yc*aux1
               fpz(i)=fpz(i)+coef*zc*aux1
            END DO
            DO k=1,diss_list(2,1)
               i=diss_list(2,k+1)
               fpx(i)=fpx(i)-coef*xc*aux2
               fpy(i)=fpy(i)-coef*yc*aux2
               fpz(i)=fpz(i)-coef*zc*aux2
            END DO
         END IF

*=======================================================================
*--- If the solute must be unfolded ------------------------------------
*=======================================================================

      ELSE
         IF(enout .LE. rsp) THEN
            enout=rsp
         ELSE
            coef=-2.0D0*alpha*(rsp-enout)/rsp
            uconf=alpha*(rsp-enout)**2
            DO k=1,diss_list(1,1)
               i=diss_list(1,k+1)
               fpx(i)=fpx(i)+coef*xc*aux1
               fpy(i)=fpy(i)+coef*yc*aux1
               fpz(i)=fpz(i)+coef*zc*aux1
            END DO
            DO k=1,diss_list(2,1)
               i=diss_list(2,k+1)
               fpx(i)=fpx(i)-coef*xc*aux2
               fpy(i)=fpy(i)-coef*yc*aux2
               fpz(i)=fpz(i)-coef*zc*aux2
            END DO
         END IF
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
