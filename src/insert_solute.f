      SUBROUTINE insert_solute(kprint,co,xp0,yp0,zp0,xpa,ypa,zpa
     &     ,xpb,ypb,zpb,index,nbtype,radius,sigma,nbtype_slv,ntap,nmol
     &     ,nato_slv)

************************************************************************
*   Time-stamp: <97/04/09 17:45:26 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 11 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ntap,nato_slv,nmol,kprint,index(*),nbtype(*),nbtype_slv(*)
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*),xpa(*),ypa(*),zpa(*),xpb(*)
     &     ,ypb(*),zpb(*),sigma(*),radius

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER m,l,i,ma,na,j,count
      REAL*8  fact,xpi,ypi,zpi,xsp,ysp,zsp,xc,yc,zc,radi,radj,radij,aux
      LOGICAL ok
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      fact=1.0D0/(2.0D0**(1.0D0/6.0D0))
      write(kprint,1)


      count=0
      DO m=1,nmol
         l=nato_slv*(m-1)
         DO i=1,nato_slv
            xpi=xpa(l+i)
            ypi=ypa(l+i)
            zpi=zpa(l+i)
            ok=.TRUE.
            radi=sigma(nbtype_slv(i))*fact
            na=ntap*(i-1)
            DO j=1,ntap
               radj=sigma(nbtype(j))*fact
               xsp=xpi-xp0(j)
               ysp=ypi-yp0(j)
               zsp=zpi-zp0(j)
               xsp=xsp-2.0D0*PBC(xsp)
               ysp=ysp-2.0D0*PBC(ysp)
               zsp=zsp-2.0D0*PBC(zsp)
               xc=co(1,1)*xsp+co(1,2)*ysp+co(1,3)*zsp
               yc=co(2,1)*xsp+co(2,2)*ysp+co(2,3)*zsp
               zc=co(3,1)*xsp+co(3,2)*ysp+co(3,3)*zsp
               xpb(na+j)=xc**2+yc**2+zc**2
               radij=(radi+radj)*radius
               radij=radij*radij
               IF(xpb(na+j) .LT. radij) THEN
                  ok=.FALSE.
                  GOTO 100
               END IF
            END DO
         END DO
100      CONTINUE
         
         aux=1.0D10
         ypb(m)=aux
         IF(ok) THEN
            count=count+1
            DO j=1,ntap*nato_slv
               IF(aux .GT. xpb(j)) aux=xpb(j)
            END DO
            ypb(m)=aux
         END IF
      END DO

*=======================================================================
*--- Order the solvent molecules according to the ascending order ------
*--- of distance from the solute atoms ---------------------------------
*=======================================================================

      CALL indexx(nmol,ypb,index)

      write(kprint,2) nmol-count,nmol
      
      nmol=count
      DO m=1,nmol
         ma=(index(m)-1)*nato_slv
         na=(m-1)*nato_slv
         DO j=1,nato_slv
            xpb(na+j)=xpa(ma+j)
            ypb(na+j)=ypa(ma+j)
            zpb(na+j)=zpa(ma+j)
         END DO
      END DO

      DO m=1,nmol*nato_slv
         xpa(m)=xpb(m)
         ypa(m)=ypb(m)
         zpa(m)=zpb(m)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1     FORMAT(//
     &     '             =============================================='
     &//   '             Inserting solute molecules            ------->'
     &//)

2     FORMAT(
     &     '             ',i6,' molecules over ',i6
     &     ,' have been removed '//)

      RETURN
      END
