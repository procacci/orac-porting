      SUBROUTINE plot_center(fstep,beta,xp0,yp0,zp0,ntap,nres,m1,prsymb
     &     ,charge)

************************************************************************
*                                                                      *
*     Print a snapshot of the solute molecules on unit kplot.          *
*                                                                      *
*     CO      :  Transform the coordinates to the orthogonalized  (I)  *
*                crystallographic frame.                               *
*                >> real*8  CO(3,3) <<                                 *
*     BETA    :  List of labels for each site of the solute.      (I)  *
*                >> character*7 BETA(NTAP) <<                          *
*     XP0     :  Coordinates of the solute molecules, packed by   (I)  *
*     YP0        site.                                                 *
*     ZP0        >> real*8  XP0(NTAP), ... <<                          *
*                                                                      *
*     NTAP    :  Number of sites on the solute molecules.         (I)  *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE.                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,m1
      INTEGER nres(m1,*)
      REAL*8  xp0(*),yp0(*),zp0(*),charge(*),fstep
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l
      REAL*8  xb,yb,zb,sunitc
      CHARACTER*5 bet,bet2
      CHARACTER*3 rsd

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kplot_center,2) fstep
      sunitc=DSQRT(unitc)
      k=0
      DO i=1,ntap
         j=nres(i,2)
         xb=xp0(i)
         yb=yp0(i)
         zb=zp0(i)
         bet(1:5)=beta(i)(1:5)
         CALL low_up(bet,5)
         IF(bet(5:5) .EQ. ' ') THEN
            bet2(1:1)=' '
            DO l=1,4
               bet2(l+1:l+1)=bet(l:l)
            END DO
         END IF
         rsd(1:3)=prsymb(j)(1:3)
         CALL low_up(rsd,3)
         WRITE(kplot_center,1)'ATOM ',i,bet2(1:5),rsd,nres(i,1),xb,yb
     &        ,zb,charge(i)*sunitc,DFLOAT(k)
      END DO
      WRITE(kplot_center,'(a)')'END  '
1     FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3,2f6.2)
2     FORMAT('REMARK   1 Configuration at time step ',f11.2,'      ',
     &'                 ')

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
