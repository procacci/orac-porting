      SUBROUTINE uscrpl(line,n)
      IMPLICIT none
      INTEGER n
      CHARACTER*1 line(n)
      CHARACTER*1 uscore,label(80)
      INTEGER i
      DATA uscore/'&'/

      DO i=1,n-1
          IF(line(i) .EQ. uscore) THEN
              GOTO 100
          END IF
      END DO
      RETURN

100   CONTINUE
      label(1)='/'
      DO i=1,n-1
          IF(line(i) .EQ. uscore) THEN
              label(i+1)=' '
          ELSE
              label(i+1)=line(i)
          END IF
      END DO
      DO i=1,n
          line(i)=label(i)
      END DO
          
      RETURN
      END
