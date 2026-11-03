      SUBROUTINE fpbond(ss_index,lbnd,lbond,xp0,yp0,zp0,pota
     &     ,potb,ubond_slt,ubond_slv,fpx,fpy,fpz,icnstr,mb,itask)

************************************************************************
*                                                                      *
*     Compute the contribution from bond stretching to solute          *
*     forces and energy.                                               *
*                                                                      *
*     CO      :  Transform the positions to orthogo alized        (I)  *
*                crystallographic frame.                               *
*                >> real*8 CO(3,3) <<                                  *
*     LBND    :  List of bends for the solute.                    (I)  *
*                >> integer LBNDG(3,N1) <<                             *
*     LBOND   :  Number of bonds.                                 (I)  *
*     N1      :  Physical column dimension of LBEND.              (I)  *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     NATO    :  Number of atoms which the macromolecule is       (I)  *
*                composed of.                                          *
*     POTA    :  List of potential parameters : bond              (I)  *
*                dostance.                                             *
*                >> real*8 POTA(N1) <<                                 *
*     POTB    :  List of potential parameters : bond              (I)  *
*                force constant.                                       *
*                >> real*8 POTB(N1) <<                                 *
*     UBOND   :  Potential energy contribution from bond.         (I)  *
*     FPX     :  List of forces acting on the macromolecule       (I)  *
*     FPY        atoms.                                                *
*     FPZ        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*                                                                      *
*---- Last update 05/17/94 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CEA, CE Saclay FRANCE 1994             *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================


#ifdef _OMP_
      use omp_integr, only:nthr1,array1_omp
#endif
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER lbond,ss_index(*)
      INTEGER lbnd(2,*),icnstr,mb,itask
      REAL*8  ubond_slt,ubond_slv
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*),
     x        fpx(mb),fpy(mb),fpz(mb)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,type,istride,itaskfp
      REAL*8  xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21,
     x        uux1,uux2,uuy1,uuy2,uuz1,uuz2,ubond(2)
      REAL*8  qforce

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         ubond(i)=0.0D0
      END DO

*----- bonded interactions: bond stretching

      DO i=1,lbond
#ifdef _OMP_
        if(array1_omp(i,itask).OR.icnstr.EQ.1) THEN
#endif
          la=lbnd(1,i)
          lb=lbnd(2,i)
          type=ss_index(la)
          x21=xp0(lb)-xp0(la)
          y21=yp0(lb)-yp0(la)
          z21=zp0(lb)-zp0(la)
          rs21=DSQRT(x21**2+y21**2+z21**2)
          qforce=-2.0D0*potb(i)*(pota(i)-rs21)/rs21
          uux1=qforce*x21
          uuy1=qforce*y21
          uuz1=qforce*z21
          fpx(la)=fpx(la)+uux1
          fpy(la)=fpy(la)+uuy1
          fpz(la)=fpz(la)+uuz1
          fpx(lb)=fpx(lb)-uux1
          fpy(lb)=fpy(lb)-uuy1
          fpz(lb)=fpz(lb)-uuz1
          ubond(type)=ubond(type)+potb(i)*(rs21-pota(i))**2
!          if(icnstr.eq.0) write(100+nthr+itask,43) i,la,lb,rs21,potb(i)
!     &         *(rs21-pota(i))**2,ubond(1),uux1,fpx(la),-uux1,fpx(lb)
!43        FORMAT(I8,2i10,f10.4,3(2x,2G12.3,2x))
#ifdef _OMP_
        ENDIF
#endif
      END DO
      ubond_slt=ubond(1)
      ubond_slv=ubond(2)
!      STOP
!      write(6,103) itask,nthr,ubond_slt,ubond_slv
!103   format(2i10," UBOND ",2f10.5)
!!$OMP BARRIER
!      if(icnstr.eq.0) STOP
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
