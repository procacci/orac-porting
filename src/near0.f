      LOGICAL FUNCTION near0(value)
      REAL*8 value,eps
      DATA eps/1.0d-6/
      near0=DABS(value).LE.eps
      RETURN
      END
