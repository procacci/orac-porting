      SUBROUTINE calc_gofr(nato_slt,nato_slv,type_slv,ss_index,atomp
     &     ,neighbor,intra,co,xp0,yp0,zp0,xpg,ypg,zpg,wca,whe,delrg,ntap
     &     ,ngrp,grppt,naa,krdf,ngrdon,rsqcut,nnlpp)

************************************************************************
*   Time-stamp: <97/08/21 11:16:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Dec  3 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 co(3,3),xp0(*),yp0(*),zp0(*),wca(*),whe(*),xpg(*),ypg(*)
     &     ,zpg(*)
      REAL*8 delrg,rsqcut
      INTEGER ntap,ngrp,naa,krdf(naa,*),grppt(2,*),ngrdon,nato_slt
     &     ,nato_slv,ss_index(*),type_slv(*),atomp(*)
      INTEGER nnlpp(*)
      LOGICAL neighbor,intra

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  delrgi,rsq,rsp,xpi,ypi,zpi,weight1,weight2,weight3,xd,yd
     &     ,zd,xc,yc,zc,xcut,xmap,ymap,zmap,xpgi,ypgi,zpgi
      INTEGER n,m,jj,i,j,noff,i1,j1,kbox,typei,typej,atom_slv,labeli
     &     ,labelj,labelij,inmax,inmin
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ngrdon=ngrdon+1
      delrgi=1.0D0/delrg
      xcut=rsqcut**2
      
      IF(neighbor) THEN
         n=0
         DO i1=1,ngrp
            m=nnlpp(1+n)
            xpgi=xpg(i1)
            ypgi=ypg(i1)
            zpgi=zpg(i1)
            IF(m .NE. 0) THEN
               typei=ss_index(grppt(1,i1))
               DO jj=1,m
                  j1=nnlpp(jj+1+n)
                  typej=ss_index(grppt(1,j1))
                  xd=xpgi-xpg(j1)
                  yd=ypgi-ypg(j1)
                  zd=zpgi-zpg(j1)
                  xmap=-2.0D0*PBC(xd)
                  ymap=-2.0D0*PBC(yd)
                  zmap=-2.0D0*PBC(zd)
                  IF(typei .EQ. 2 .AND. typej .EQ. 2) THEN
                     DO i=grppt(1,i1),grppt(2,i1)
                        xpi=xp0(i)
                        ypi=yp0(i)
                        zpi=zp0(i)
                        atom_slv=MOD(i-nato_slt-1,nato_slv)+1
                        labeli=type_slv(atom_slv)
                        IF(atomp(i) .NE. atomp(grppt(1,j1)) .OR. intra)
     &                       THEN
                           DO j=grppt(1,j1),grppt(2,j1)
                              atom_slv=MOD(j-nato_slt-1,nato_slv)+1
                              labelj=type_slv(atom_slv)
                              inmax=max0(labeli,labelj)
                              inmin=min0(labeli,labelj)
                              labelij=inmax*(inmax-1)/2+inmin
                              xd=xpi-xp0(j)
                              yd=ypi-yp0(j)
                              zd=zpi-zp0(j)
                              xd=xd+xmap
                              yd=yd+ymap
                              zd=zd+zmap
                              xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                              yc=           co(2,2)*yd+co(2,3)*zd
                              zc=                      co(3,3)*zd
                              rsq=xc*xc+yc*yc+zc*zc
                              rsp=DSQRT(rsq)
                              kbox=MIN0(IDINT(delrgi*rsp+0.5D0),maxint)
                              krdf(kbox,3+labelij)=krdf(kbox,3+labelij)
     &                             +1
                           END DO
                        END IF
                     END DO
                  ELSE IF(typei .EQ. 1 .AND. typej .EQ. 1) THEN
                     DO i=grppt(1,i1),grppt(2,i1)
                        xpi=xp0(i)
                        ypi=yp0(i)
                        zpi=zp0(i)
                        DO j=grppt(1,j1),grppt(2,j1)
                           weight1=wca(i)*wca(j)
                           weight2=whe(i)*whe(j)
                           weight3=1.0D0
                           xd=xpi-xp0(j)
                           yd=ypi-yp0(j)
                           zd=zpi-zp0(j)
                           xd=xd+xmap
                           yd=yd+ymap
                           zd=zd+zmap
                           xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                           yc=           co(2,2)*yd+co(2,3)*zd
                           zc=                      co(3,3)*zd
                           rsq=xc*xc+yc*yc+zc*zc
                           rsp=DSQRT(rsq)
                           kbox=MIN0(IDINT(delrgi*rsp+0.5D0),maxint)
                           krdf(kbox,1)=krdf(kbox,1)+IDINT(weight1)
                           krdf(kbox,2)=krdf(kbox,2)+IDINT(weight2)
                           krdf(kbox,3)=krdf(kbox,3)+IDINT(weight3)
                        END DO
                     END DO
                  END IF
               END DO
            END IF
            noff=nnlpp(1+n)+1
            n=n+noff
         END DO
      ELSE
         DO i=1,ntap-1
            xpi=xp0(i)
            ypi=yp0(i)
            zpi=zp0(i)
            typei=ss_index(i)
            DO j=i+1,ntap
               typej=ss_index(j)
               IF(atomp(i) .NE. atomp(j) .OR. intra) THEN
                  IF(typei .EQ. 2 .AND. typej .EQ. 2) THEN
                     atom_slv=MOD(i-nato_slt-1,nato_slv)+1
                     labeli=type_slv(atom_slv)
                     atom_slv=MOD(j-nato_slt-1,nato_slv)+1
                     labelj=type_slv(atom_slv)
                     inmax=max0(labeli,labelj)
                     inmin=min0(labeli,labelj)
                     labelij=inmax*(inmax-1)/2+inmin
                     xd=xpi-xp0(j)
                     yd=ypi-yp0(j)
                     zd=zpi-zp0(j)
                     xd=xd-2.0D0*PBC(xd)
                     yd=yd-2.0D0*PBC(yd)
                     zd=zd-2.0D0*PBC(zd)
                     xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                     yc=           co(2,2)*yd+co(2,3)*zd
                     zc=                      co(3,3)*zd
                     rsq=xc*xc+yc*yc+zc*zc
                     IF(xcut .GT. rsq) THEN
                        rsp=DSQRT(rsq)
                        kbox=MIN0(IDINT(delrgi*rsp+0.5D0),maxint)
                        krdf(kbox,3+labelij)=krdf(kbox,3+labelij)+1
                     END IF
                  END IF
               ELSE IF(typei .EQ. 1 .AND. typej .EQ. 1) THEN
                  weight1=wca(i)*wca(j)
                  weight2=whe(i)*whe(j)
                  weight3=1.0D0
                  xd=xpi-xp0(j)
                  yd=ypi-yp0(j)
                  zd=zpi-zp0(j)
                  xd=xd-2.0D0*PBC(xd)
                  yd=yd-2.0D0*PBC(yd)
                  zd=zd-2.0D0*PBC(zd)
                  xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                  yc=co(2,1)*xd+co(2,2)*yd+co(2,3)*zd
                  zc=co(3,1)*xd+co(3,2)*yd+co(3,3)*zd
                  rsq=xc*xc+yc*yc+zc*zc
                  IF(xcut .GT. rsq) THEN
                     rsp=DSQRT(rsq)
                     kbox=MIN0(IDINT(delrgi*rsp+0.5D0),maxint)
                     krdf(kbox,1)=krdf(kbox,1)+IDINT(weight1)
                     krdf(kbox,2)=krdf(kbox,2)+IDINT(weight2)
                     krdf(kbox,3)=krdf(kbox,3)+IDINT(weight3)
                  END IF
               END IF
            END DO
         END DO
      END IF

      IF(intra) THEN
         DO i1=1,ngrp
            typei=ss_index(grppt(1,i1))
            IF(typei .EQ. 2) THEN
               DO i=grppt(1,i1),grppt(2,i1)-1
                  xpi=xp0(i)
                  ypi=yp0(i)
                  zpi=zp0(i)
                  atom_slv=MOD(i-nato_slt-1,nato_slv)+1
                  labeli=type_slv(atom_slv)
                  DO j=i+1,grppt(2,i1)
                     atom_slv=MOD(j-nato_slt-1,nato_slv)+1
                     labelj=type_slv(atom_slv)
                     inmax=max0(labeli,labelj)
                     inmin=min0(labeli,labelj)
                     labelij=inmax*(inmax-1)/2+inmin
                     xd=xpi-xp0(j)
                     yd=ypi-yp0(j)
                     zd=zpi-zp0(j)
                     xd=xd+xmap
                     yd=yd+ymap
                     zd=zd+zmap
                     xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                     yc=           co(2,2)*yd+co(2,3)*zd
                     zc=                      co(3,3)*zd
                     rsq=xc*xc+yc*yc+zc*zc
                     rsp=DSQRT(rsq)
                     kbox=MIN0(IDINT(delrgi*rsp+0.5D0),maxint)
                     krdf(kbox,3+labelij)=krdf(kbox,3+labelij)
     &                    +1
                  END DO
               END DO
            END IF
         END DO
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
