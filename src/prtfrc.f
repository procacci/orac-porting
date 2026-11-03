      SUBROUTINE prtfrc(kprint,ngrp,grppt,nres,M1,prsymb,beta,xp0,yp0
     &     ,zp0,fpx,fpy,fpz)
      IMPLICIT none
      INTEGER kprint,ngrp,m1,grppt(2,*),nres(m1,*)
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
      CHARACTER*8 prsymb(*)
      CHARACTER*7 beta(*)
     
      INTEGER i,i1

      DO i=1,ngrp
         DO i1=grppt(1,i),grppt(2,i)
            WRITE(kprint
     &           ,'(5x,3(f10.5,1x),3(e12.5,1x),2i6,1x,a7,1x,a8,i6)'
     &           )xp0(i1),yp0(i1),zp0(i1),fpx(i1),fpy(i1),fpz(i1),i1,i
     &           ,beta(i1),prsymb(nres(i1,2)),nres(i1,1)
         END DO
      END DO
      RETURN
      END
