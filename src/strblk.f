      INTEGER FUNCTION strblk(line,n)
      IMPLICIT none
      INTEGER n
      CHARACTER*1 line(n)
      INTEGER i

      strblk=0
      DO i=1,n
          IF(line(i) .EQ. ' ' .OR. line(i) .EQ. ' ') THEN
              strblk=i
              RETURN
          END IF
      END DO
      RETURN
      END
