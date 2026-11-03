      RECURSIVE SUBROUTINE recpr(index,concta,m1,mask)
      IMPLICIT none
      INTEGER index,m1
      INTEGER concta(m1,*)
      LOGICAL mask(*)

      INTEGER i,coordi,ia

      mask(index)=.FALSE.
      coordi=concta(index,1)
      DO i=2,coordi+1
          ia=concta(index,i)
          IF(mask(ia)) THEN
              CALL recpr(ia,concta,m1,mask)
          END IF
      END DO
      RETURN
      END
