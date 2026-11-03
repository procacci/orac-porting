      SUBROUTINE fpbend(ss_index,lbndg,lbend,xp0,yp0,zp0,pota,potb
     &     ,potc,potd,pote,potf,potg,ubend_slt,ubend_slv,fpx,fpy,fpz,mb
     &     ,itask)
************************************************************************
*                                                                      *
*     Compute the contribution from bond angles interaction            *
*     to macromolecular forces and energy.                             *
*                                                                      *
*     CO      :  Transform the positions to orthogo alized        (I)  *
*                crystallographic frame.                               *
*                >> real*8 CO(3,3) <<                                  *
*     LBNDG   :  List of bendings for the macromolecule.          (I)  *
*                >> integer LBNDG(3,N1) <<                             *
*     LBEND   :  Number of bendings.                              (I)  *
*     N1      :  Physical column dimension of LBEND.              (I)  *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     NATO    :  Number of atoms which the macromolecule is       (I)  *
*                composed of.                                          *
*     POTA    :  List of potential parameters : bending           (I)  *

*                theta angles.                                         *
*                >> real*8 POTA(N1) <<                                 *
*     POTB    :  List of potential parameters : bending           (I)  *
*                force constant.                                       *
*                >> real*8 POTB(N1) <<                                 *
*     UBEND   :  Potential energy contribution from bendings.     (I)  *
*     FPX     :  List of forces acting on the macromolecule       (I)  *
*     FPY        atoms.                                                *
*     FPZ        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*                                                                      *
*---- Last update 04/24/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

#ifdef _OMP_
      use omp_integr, only:nthr1,array1_omp
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------



      INTEGER lbend,ss_index(*),mb
      INTEGER lbndg(3,*),itask
      REAL*8  ubend_slt,ubend_slv
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*),potc(*),potd(*),pote(
     &     *),potf(*),potg(*),fpx(mb),fpy(mb),fpz(mb)
*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,lc,type
      REAL*8  xr1,xr2,xr3,yr1,yr2,yr3,zr1,zr2,zr3,x12,x32,y12,y32,z12
     &     ,z32,rs12,rs32,uux1,uux2,uux3,uuy1,uuy2,uuy3,uuz1,uuz2,uuz3
     &     ,xr31,yr31,zr31,rs31,rsp31,x31,y31,z31,rsp12,rsp32,k12,r1,r2
      REAL*8  dcc2,cb,sb,bb,qforce,pforce,pi,ubend(2)
      logical coupling 

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         ubend(i)=0.0D0
      END DO
      pi=4.0D0*DATAN(1.0D0)

*----- bonded interactions: bending

      DO i=1,lbend
         if(abs(potb(i)).lt.1.d-12) go to 1002
#ifdef _OMP_
         if(array1_omp(i,itask)) THEN
#endif
         la=lbndg(1,i)
         lb=lbndg(2,i)
         lc=lbndg(3,i)
         type=ss_index(la)
         xr1=xp0(la)
         yr1=yp0(la)
         zr1=zp0(la)
         xr2=xp0(lb)
         yr2=yp0(lb)
         zr2=zp0(lb)
         xr3=xp0(lc)
         yr3=yp0(lc)
         zr3=zp0(lc)
         x12=xr1-xr2
         y12=yr1-yr2
         z12=zr1-zr2
         x32=xr3-xr2
         y32=yr3-yr2
         z32=zr3-zr2
         rs12=x12**2+y12**2+z12**2
         rs32=x32**2+y32**2+z32**2
         dcc2=DSQRT(rs12*rs32)
         cb=(x12*x32+y12*y32+z12*z32)/dcc2
         if(dabs(cb).gt.1.0) cb=sign(1.d0,cb) 
         sb=DSQRT(1.0d0-cb**2)
         bb=DACOS(cb)
c--      switch to singularity free potential for 
c--      linear bending V = 2*K*(cos(theta) + 1)  
         if(abs(pota(i)-pi).lt.0.01.or.sb.lt.1D-6) THEN 
           qforce=-2.0d0*potb(i)
c           write(6,77) i,bb*180.d0/pi,qforce,2.0d0*potb(i)*(bb-pota(i))
c     &          /sb
         else
           qforce=2.0d0*potb(i)*(bb-pota(i))/sb
         end if
         uux1=x32/dcc2-cb*x12/rs12
         uux3=x12/dcc2-cb*x32/rs32
         uux2=-uux1-uux3
         uuy1=y32/dcc2-cb*y12/rs12
         uuy3=y12/dcc2-cb*y32/rs32
         uuy2=-uuy1-uuy3
         uuz1=z32/dcc2-cb*z12/rs12
         uuz3=z12/dcc2-cb*z32/rs32
         uuz2=-uuz1-uuz3
         fpx(la)=fpx(la)+qforce*uux1
         fpx(lc)=fpx(lc)+qforce*uux3
         fpx(lb)=fpx(lb)+qforce*uux2
         fpy(la)=fpy(la)+qforce*uuy1
         fpy(lc)=fpy(lc)+qforce*uuy3
         fpy(lb)=fpy(lb)+qforce*uuy2
         fpz(la)=fpz(la)+qforce*uuz1
         fpz(lc)=fpz(lc)+qforce*uuz3
         fpz(lb)=fpz(lb)+qforce*uuz2
c--      switch to singularity free potential for 
c--      linear bending V = 2*K(cos(theta) + 1)  
         if(abs(pota(i)-pi).lt.0.01) THEN 
           ubend(type)=ubend(type)+2.d0*potb(i)*(cb + 1.d0) 
c           write(6,77) i,bb*180.d0/pi,qforce,2.0d0*potb(i)*(bb-pota(i))/sb
c     &          ,2.d0*potb(i)*(cb + 1.d0),potb(i)*(bb-pota(i))**2
c77         FORMAT(i10,f12.6,4e15.5)
         else
           ubend(type)=ubend(type)+potb(i)*(bb-pota(i))**2
         end if
*=======================================================================
*---  Compute Coupling terms
*=======================================================================

c        the Urey-Bradley term 
c        accept syntax   kub  rub 0.0 0.0 0.0 or 
c        accept syntax   kub  rub k12 r1 r2
         if(abs(potd(i)).gt.1.d-5.and.(abs(pote(i)).lt.1.d-5.or
     &        .abs(potf(i)).gt.1.d-5)) THEN 
           x31=xr3-xr1
           y31=yr3-yr1
           z31=zr3-zr1
           rs31=x31**2+y31**2+z31**2
           rsp31=DSQRT(rs31)
           
           qforce=-2.0D0*potd(i)*(potc(i)-rsp31)
           uux1=x31/rsp31
           uuy1=y31/rsp31
           uuz1=z31/rsp31
           uux2=-uux1
           uuy2=-uuy1
           uuz2=-uuz1
           fpx(la)=fpx(la)+qforce*uux1
           fpy(la)=fpy(la)+qforce*uuy1
           fpz(la)=fpz(la)+qforce*uuz1
           fpx(lc)=fpx(lc)+qforce*uux2
           fpy(lc)=fpy(lc)+qforce*uuy2
           fpz(lc)=fpz(lc)+qforce*uuz2
           ubend(type)=ubend(type)+potd(i)*(rsp31-potc(i))**2 
         END IF
c        find out whether coupling for bending is present 
         coupling = .false. 
         IF((abs(potd(i)).gt.1.d-5.and.abs(pote(i)).gt.1.d-5.and
     &        .abs(potf(i)).lt.1.d-5)) THEN 
           coupling = .true. 
           k12 = potd(i) 
           r1 = potc(i) 
           r2 = pote(i) 
         ENDIF 
         IF((abs(pote(i)).gt.1.d-5.and.abs(potf(i)).gt.1.d-5.and
     &        .abs(potg(i)).GT.1.d-5)) THEN 
           coupling = .true. 
           k12 = pote(i) 
           r1 = potf(i) 
           r2 = potg(i) 
         ENDIF 
         IF(coupling) THEN
c--      off diagonal coupling for stretchings 
           rsp12 = dsqrt(rs12) 
           rsp32 = dsqrt(rs32) 
           uux1=x12/rsp12
           uuy1=y12/rsp12
           uuz1=z12/rsp12
           uux2=x32/rsp32
           uuy2=y32/rsp32
           uuz2=z32/rsp32
           qforce = k12*(rsp12-r1) 
           pforce = k12*(rsp32-r2) 
           fpx(la)=fpx(la)-pforce*uux1
           fpy(la)=fpy(la)-pforce*uuy1
           fpz(la)=fpz(la)-pforce*uuz1
           fpx(lb)=fpx(lb) + pforce*uux1 + qforce*uux2
           fpy(lb)=fpy(lb) + pforce*uuy1 + qforce*uuy2
           fpz(lb)=fpz(lb) + pforce*uuz1 + qforce*uuz2
           fpx(lc)=fpx(lc)- qforce*uux2
           fpy(lc)=fpy(lc)- qforce*uuy2
           fpz(lc)=fpz(lc)- qforce*uuz2
           ubend(type)=ubend(type)+k12*(rsp12-r1)
     &          *(rsp32-r2) 
         END IF
#ifdef _OMP_
         ENDIF
#endif
1002     CONTINUE
      END DO

      ubend_slt=ubend(1)
      ubend_slv=ubend(2)
!      write(6,103) itask,nthr,ubend_slt,ubend_slv
!103   format(2i10," UBEND ",2f10.5)
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
