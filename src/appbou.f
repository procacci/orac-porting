      SUBROUTINE appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)

************************************************************************
*                                                                      *
*     Recalculate positions of the groups                              *
*                                                                      *
*     XP0     :  Protein coordinates in simulation box units.   (I/O)  *
*     YP0        >> REAL*8 XP0(*) <<                                   *
*     ZP0                                                              *
*     XPG     :  Coordinates of the groups in simulation box    (I/O)  *
*     YPG        units.                                                *
*     ZPG        >> REAL*8 XPG(*) <<                                   *
*     PMASS   :  List of atomic masses.                           (I)  *
*                >> REAL*8 PMASS(*) <<                                 *
*     NGRP    :  Number of groups in the protein                  (I)  *
*     GRPPT   :  Pointer to the first an last atom of a group     (I)  *
*                >> INTEGER GRPPT(2,*) <<                              *
*     RGGRUP  :  If SECTION is .TRUE. this is the number of       (I)  *
*                groups which are contained in the primary             *
*                region                                                *
*     RGGMAP  :  If SECTION is .TRUE. this is the list of atom    (I)  *
*                numbers contained in the primary region               *
*                >> INTEGER RGGMAP(*) <<                               *
*     SCTION  :  Flag which indicates if a primary and secondary  (I)  *
*                simulation regions have been defined                  *
*                >> LOGICAL SCTION <<                                  *
*                                                                      *
*---- Last update  25/1092 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992               *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use omp_integr, only:omp_timing,time_in,tcpu_in
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ngrp,grppt(2,*)
      REAL*8  pmass(*),xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  sumx,sumy,sumz
      INTEGER n,m
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal

*==================== EXECUTABLE STATEMENTS ============================

      if(omp_timing)  call timer(treal,tcpu1,tela1)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(sumx,sumy,sumz)
!$OMP& SHARED(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp)
      DO n=1,ngrp
         sumx=0.0D0
         sumy=0.0D0
         sumz=0.0D0
         DO m=grppt(1,n),grppt(2,n)
            sumx=sumx+pmass(m)*xp0(m)
            sumy=sumy+pmass(m)*yp0(m)
            sumz=sumz+pmass(m)*zp0(m)
         END DO
         xpg(n)=sumx
         ypg(n)=sumy
         zpg(n)=sumz
      END DO
!$OMP END PARALLEL DO
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_in=time_in+tela2-tela1
        tcpu_in=tcpu_in+tcpu2-tcpu1
      END IF


*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
