      SUBROUTINE update(co,xpg,ypg,zpg,xp0,yp0,zp0,rspcut
     &     ,rspoff,ngrp,nnlpp,o1)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 12 1997 -                                     *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ngrp,o1
      INTEGER nnlpp(*)
      REAL*8  co(3,3),xpg(*),ypg(*),zpg(*),xp0(*),yp0(*),zp0(*),rspcut
     &     ,rspoff

*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nb
      PARAMETER (nb=m11)
      REAL*8  xmap(nb),ymap(nb),zmap(nb)
      INTEGER mmap(nb)
      LOGICAL maplg(2*nb)
      COMMON / rag1 / xmap,ymap,zmap,mmap,maplg

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER i,j,m,n,ia,ira,id,ird,map,ncount,l1,l2,
     x        l3,l4
      REAL*8  radius,x32,y32,z32,xc32,yc32,zc32,rsq32,rsp32
      REAL*8  xpi,ypi,zpi,xsa,ysa,zsa,xca,yca,zca
      INTEGER  npp
      INCLUDE 'pbc.h'

*==================== EXECUTABLE STATEMENTS ============================

      IF(nb .LT. ngrp) THEN
          WRITE(6,*) nb,m11
          WRITE(6,*) ngrp
          errmsg=' IN UPDATE : Dimensions of the work arrays are '
     &/ /'insufficient.'
          CALL xerror(errmsg,80,1,2)
          STOP
      END IF


*=======================================================================
*----- THE SIMULATION INVOLVES SOLUTE ATOMS ----------------------------
*=======================================================================

      radius=(rspcut+rspoff)**2
      ncount=0
      n=0
      DO i=1,ngrp
         xpi=xpg(i)
         ypi=ypg(i)
         zpi=zpg(i)

*=======================================================================
*----- Fill up the position maps accordingly ---------------------------
*=======================================================================

         DO j=i+1,ngrp
            xsa=xpi-xpg(j)-2.0D0*PBC(xpi-xpg(j))
            ysa=ypi-ypg(j)-2.0D0*PBC(ypi-ypg(j))
            zsa=zpi-zpg(j)-2.0D0*PBC(zpi-zpg(j))
            xca=co(1,1)*xsa+co(1,2)*ysa+co(1,3)*zsa
            yca=co(2,1)*xsa+co(2,2)*ysa+co(2,3)*zsa
            zca=co(3,1)*xsa+co(3,2)*ysa+co(3,3)*zsa
            xmap(j)=xca*xca+yca*yca+zca*zca
         END DO

*=======================================================================
*----- Take the atoms which are inside sphere of radius RADIUS ---------
*=======================================================================
         
         m=0
         DO j=i+1,ngrp
            IF(xmap(j) .LE. radius) THEN
               m=m+1
               nnlpp(m+1+ncount)=j
            END IF
         END DO
         
         nnlpp(ncount+1)=m
         ncount=ncount+nnlpp(ncount+1)+1
         IF(o1.LT.ncount) THEN
            errmsg=' IN UPDATE : Dimensions of NNLPP '/ /
     x           'are insufficient. ABORT !! '
            WRITE(kprint,'(a)') errmsg
            WRITE(kprint,'('' O1 ='',i7,'' NCOUNT = '',i7)') o1,ncount
            CALL xerror(errmsg,80,1,2)
         END IF
      END DO

      npp=ncount

      WRITE(kprint,10000) npp

*================= END OF EXECUTABLE STATEMENTS ========================

10000 FORMAT('Neighbor Lists Dimensions     *neighbor(',i7,')* '/)

      RETURN
      END
