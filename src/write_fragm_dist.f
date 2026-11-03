      SUBROUTINE write_fragm_dist(fstep,kfrag,fragm,nfrag,xpa,ypa,zpa,co
     &     )

************************************************************************
*   Time-stamp: <97/07/09 15:08:36 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  9 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER fragm(2,*),nfrag,kfrag
      REAL*8  xpa(*),ypa(*),zpa(*),co(3,3),fstep

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,i1,j,j1
      REAL*8  xpfi,ypfi,zpfi,xpfj,ypfj,zpfj,xd,yd,zd,xc,yc,zc,aux,rsq
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO i=1,nfrag
         xpfi=0.0D0
         ypfi=0.0D0
         zpfi=0.0D0
         DO i1=fragm(1,i),fragm(2,i)
            xpfi=xpfi+xpa(i1)
            ypfi=ypfi+ypa(i1)
            zpfi=zpfi+zpa(i1)
         END DO
         aux=1.0D0/DFLOAT(fragm(2,i)-fragm(1,i)+1)
         xpfi=xpfi*aux
         ypfi=ypfi*aux
         zpfi=zpfi*aux
         DO j=i+1,nfrag
            xpfj=0.0D0
            ypfj=0.0D0
            zpfj=0.0D0
            DO j1=fragm(1,j),fragm(2,j)
               xpfj=xpfj+xpa(j1)
               ypfj=ypfj+ypa(j1)
               zpfj=zpfj+zpa(j1)
            END DO
            aux=1.0D0/DFLOAT(fragm(2,j)-fragm(1,j)+1)
            xpfj=xpfj*aux
            ypfj=ypfj*aux
            zpfj=zpfj*aux
            xd=xpfi-xpfj
            yd=ypfi-ypfj
            zd=zpfi-zpfj
            xd=xd-2.0D0*PBC(xd)
            yd=yd-2.0D0*PBC(yd)
            zd=zd-2.0D0*PBC(zd)
            xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
            yc=co(2,1)*xd+co(2,2)*yd+co(2,3)*zd
            zc=co(3,1)*xd+co(3,2)*yd+co(3,3)*zd
            rsq=xc*xc+yc*yc+zc*zc
            WRITE(kfrag,100) 'T',fstep,DSQRT(rsq),i,j
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(a1,f11.2,3x,f12.5,3x,2i4)
      RETURN
      END
