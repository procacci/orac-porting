      SUBROUTINE plotd(fstep,kout,beta,xp_ini,yp_ini,zp_ini,xp0,yp0,zp0
     &     ,ntap,nres,m1,prsymb)

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

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,m1,kout
      INTEGER nres(m1,*)
      REAL*8  xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*),zp_ini(*),fstep
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l,m,ndx
      REAL*8  xb,yb,zb,charge,dx,dy,dz,dr,xd,yd,zd,sum
      CHARACTER*5 bet,bet2
      CHARACTER*3 rsd

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kout,2) fstep
      charge=0.0D0
      k=0
      DO i=1,ntap
         j=nres(i,2)
         xb=xp0(i)
         yb=yp0(i)
         zb=zp0(i)
         sum=DABS(xb)+DABS(yb)+DABS(zb)
         IF(sum .GT. 1.0D-3) THEN
            IF(xp_ini(i) .GT. 1.0D4) THEN
               dr=0.0D0
            ELSE
               dx=xb-xp_ini(i)
               dy=yb-yp_ini(i)
               dz=zb-zp_ini(i)
               dr=DSQRT(dx**2+dy**2+dz**2)
            END IF
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
            WRITE(kout,1)'ATOM ',i,bet2(1:5),rsd,nres(i,1),xb,yb,zb
     &           ,dr,DFLOAT(k)
         END IF
      END DO

1     FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3,f8.4,f4.1)
2     FORMAT('REMARK   1 Configuration at time step ',f11.2,'      ',
     &'                 ')

*================= END OF EXECUTABLE STATEMENTS ========================
      
      RETURN
      END
