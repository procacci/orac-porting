      SUBROUTINE repl_coord(xa,ya,za,na,natoa,cg,tg,gnmol,icl,icm,icn
     &     ,nmol,nts,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/12 21:49:21 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 13 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  xa(*),ya(*),za(*),cg(3,3,*),tg(3,*)
      INTEGER na,icl,icm,icn,nmol,natoa,nts,gnmol,iret
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,l,n,m,count,countm
      REAL*8  dl,dm,dn,xd,yd,zd

*----------------------- EXECUTABLE STATEMENTS ------------------------*


!=======================================================================
!--------- Set up the elementary cell ----------------------------------
!=======================================================================

      countm=1
      count=natoa
      DO n=2,gnmol
         DO m=1,natoa
            l=natoa*(n-1)+m
            xd=cg(1,1,n)*xa(m)+cg(2,1,n)*ya(m)+cg(3,1,n)*za(m)
            yd=cg(1,2,n)*xa(m)+cg(2,2,n)*ya(m)+cg(3,2,n)*za(m)
            zd=cg(1,3,n)*xa(m)+cg(2,3,n)*ya(m)+cg(3,3,n)*za(m)
            xa(l)=xd+boxl*tg(1,n)/DFLOAT(icl)
            ya(l)=yd+boxl*tg(2,n)/DFLOAT(icm)
            za(l)=zd+boxl*tg(3,n)/DFLOAT(icn)
         END DO
         count=count+natoa
         countm=countm+1
      END DO
      nmol=0
      nts=0
      DO l=1,icl
         dl=DFLOAT(l-1)/DFLOAT(icl)
         DO m=1,icm
            dm=DFLOAT(m-1)/DFLOAT(icm)
            DO n=1,icn
               dn=DFLOAT(n-1)/DFLOAT(icn)
               DO j=1,count
                  nts=nts+1
                  i=nts
                  IF(i .GT. na) THEN
                     errmsg=
     &'In REPL_COORD: Coordinate dimensions insufficient. Abort.'
                     iret=1
                     WRITE(6,*) na,i
                     RETURN
                  END IF
                  xa(i)=xa(j) + boxl*dl
                  ya(i)=ya(j) + boxl*dm
                  za(i)=za(j) + boxl*dn
               END DO
               nmol=nmol+countm
            END DO
         END DO
      END DO
      nts=i

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
