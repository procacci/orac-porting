      SUBROUTINE pntres(res,pointer,atomg,grppt,nunit,resg,ngrp,nato)
      IMPLICIT none

      INTEGER res(*),resg(*),pointer(2,*),grppt(2,*),atomg(*),nunit,nato
     &     ,ngrp

      INTEGER n,m,i,ia

      n=1
      pointer(1,1)=1
      ia=res(1)
      resg(1)=ia
      DO m=1,ngrp
         DO i=grppt(1,m),grppt(2,m)
            IF(res(i) .NE. ia) THEN
               ia=res(i)
               pointer(2,n)=m-1
               n=n+1
               pointer(1,n)=m
            END IF
            atomg(i)=m
         END DO
         resg(m)=ia
      END DO
      pointer(2,n)=ngrp
      RETURN
      END
