      SUBROUTINE chlab(subcmd,res,label,labelb,nlbl,nlblb,ndimb,
     x                 iret,errmsg)
      IMPLICIT none
      INTEGER nlbl,nlblb,ndimb,iret
      CHARACTER*7 label(*),labelb(ndimb,*)
      CHARACTER*8 res,subcmd
      CHARACTER*80 errmsg

      INTEGER n,m,l
      LOGICAL ok

      iret=0
      errmsg=' '
      DO n=1,nlblb
          DO m=1,ndimb
              IF(labelb(m,n)(1:1) .NE. '-' .AND. labelb(m,n)(1:1)
     x           .NE. '+' .AND. labelb(m,n)(1:1) .NE. '*') THEN
                  ok=.FALSE.
                  DO l=1,nlbl
                      IF(label(l) .EQ. labelb(m,n)) THEN
                          ok=.TRUE.
                      END IF
                  END DO
                  IF(.NOT. ok) THEN
                      iret=1
                      errmsg=
     x'In Chlab for subcmd '//subcmd//': Label '//labelb(m,n)//
     x'in res '//res(1:5)//' not found in atom list. Abort.'
                      RETURN
                  END IF
              END IF
          END DO
      END DO
      RETURN
      END
