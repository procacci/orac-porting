      SUBROUTINE redmss(ngrp,grppt,mass,pmass)
      IMPLICIT none
      INTEGER ngrp,grppt(2,*)
      REAL*8  mass(*),pmass(*)
      
      INTEGER n,m
      REAL*8  sum

      DO n=1,ngrp
          sum=0.0D0
          DO m=grppt(1,n),grppt(2,n)
              sum=sum+mass(m)
          END DO
          DO m=grppt(1,n),grppt(2,n)
              pmass(m)=mass(m)/sum
          END DO
      END DO
      RETURN
      END
