      SUBROUTINE generate_slv(co,oc,xpa,ypa,zpa,xpb,ypb,zpb,mass,nato_in
     &     ,nmol,icl,icm,icn,nform,rmol,boxl,slv_randomize,iret,errmsg)

************************************************************************
*   Time-stamp: <97/02/28 10:53:08 marchi>                             *
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

      INTEGER nato_in,icl,icm,icn,nform,nmol,iret
      REAL*8  xpa(*),ypa(*),zpa(*),xpb(*),ypb(*),zpb(*),rmol(3,*),mass(
     &     *),boxl,oc(3,3),co(3,3)
      CHARACTER*80 errmsg
      LOGICAL slv_randomize

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,l,m,n,k,j
      REAL*8  shiftl,shiftm,shiftn,dl,dm,dn,shifti,sum,xcm,ycm,zcm,fact
     &     ,x0i,y0i,z0i,d(3,3),xknock,yknock,zknock,bumpx,bumpy,bumpz
     &     ,ranf
      external ranf

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      sum=0.0D0
      xcm=0.0D0
      ycm=0.0D0
      zcm=0.0D0
      
      DO i=1,nato_in
         fact=mass(i)
         xcm=xcm+fact*xpa(i)
         ycm=ycm+fact*ypa(i)
         zcm=zcm+fact*zpa(i)
         sum=sum+fact
      END DO
      xcm=xcm/sum
      ycm=ycm/sum
      zcm=zcm/sum
      DO i=1,nato_in
         xpa(i)=xpa(i)-xcm
         ypa(i)=ypa(i)-ycm
         zpa(i)=zpa(i)-zcm
      END DO
*-----------------------------------------------------------------------      
*---- Generate coordinates of the center of mass over the chosen grid --
*-----------------------------------------------------------------------      
      xknock=(oc(1,1)+oc(1,2)+oc(1,3))*0.3D0
      yknock=(oc(1,1)+oc(1,2)+oc(1,3))*0.3D0
      zknock=(oc(1,1)+oc(1,2)+oc(1,3))*0.3D0

      shiftl=1.0d0/DFLOAT(icl)
      shiftm=1.0d0/DFLOAT(icm)
      shiftn=1.0d0/DFLOAT(icn)
      shifti=boxl

      i=0
      DO l=1,icl
         dl=DFLOAT(l)-0.5D0
         DO m=1,icm
            dm=DFLOAT(m)-0.5D0
            DO n=1,icn
               dn=DFLOAT(n)-0.5D0
               DO j=1,nform
                  i=i+1
                  x0i=shifti*(shiftl*(rmol(1,j) + dl)-0.5D0)
                  y0i=shifti*(shiftm*(rmol(2,j) + dm)-0.5D0)
                  z0i=shifti*(shiftn*(rmol(3,j) + dn)-0.5D0)
                  bumpx=2.0D0*(ranf()-1.0D0)*xknock
                  bumpy=2.0D0*(ranf()-1.0D0)*yknock
                  bumpz=2.0D0*(ranf()-1.0D0)*zknock
                  IF(slv_randomize) THEN
                     xpb(i)=x0i+bumpx
                     ypb(i)=y0i+bumpy
                     zpb(i)=z0i+bumpz
                  ELSE
                     xpb(i)=x0i
                     ypb(i)=y0i
                     zpb(i)=z0i
                  END IF
               END DO
            END DO
         END DO
      END DO
      
*-----------------------------------------------------------------------      
*---- Compute first solvent molecule center of mass --------------------
*-----------------------------------------------------------------------      

      IF(i.NE.nmol) THEN
         iret=1
         errmsg=
     &        'While generating solvent: Number of generated molecule'
     &      //' does not match expected number.'
         RETURN
      END IF
      DO i=1,3
         DO j=1,3
            IF(i .NE. j) THEN
               d(i,j)=0.0D0
            ELSE
               d(i,j)=1.0D0
            END IF
         END DO
      END DO
      

      DO m=1,nmol
         IF(slv_randomize) THEN
            CALL generate_rnd_rot(d)
         END IF
         n=nato_in*m
         DO k=1,nato_in
            xcm=d(1,1)*xpa(k)+d(1,2)*ypa(k)+d(1,3)*zpa(k)
            ycm=d(2,1)*xpa(k)+d(2,2)*ypa(k)+d(2,3)*zpa(k)
            zcm=d(3,1)*xpa(k)+d(3,2)*ypa(k)+d(3,3)*zpa(k)
            xpa(n+k)=xpb(m)+oc(1,1)*xcm+oc(1,2)*ycm+oc(1,3)*zcm
            ypa(n+k)=ypb(m)+oc(2,1)*xcm+oc(2,2)*ycm+oc(2,3)*zcm
            zpa(n+k)=zpb(m)+oc(3,1)*xcm+oc(3,2)*ycm+oc(3,3)*zcm
         END DO
      END DO

      DO m=1,nmol
         n=nato_in*m
         DO k=1,nato_in
            xpa(n-nato_in+k)=xpa(n+k)
            ypa(n-nato_in+k)=ypa(n+k)
            zpa(n-nato_in+k)=zpa(n+k)
         END DO
      END DO
      n=nmol*nato_in
      DO m=1,n
         xcm=co(1,1)*xpa(m)+co(1,2)*ypa(m)+co(1,3)*zpa(m)
         ycm=co(2,1)*xpa(m)+co(2,2)*ypa(m)+co(2,3)*zpa(m)
         zcm=co(3,1)*xpa(m)+co(3,2)*ypa(m)+co(3,3)*zpa(m)
         xpa(m)=xcm
         ypa(m)=ycm
         zpa(m)=zcm
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
