      SUBROUTINE xfit(xyz,xyz0,xyzfit,qt,w1,w2,ntap,error)
      IMPLICIT none
      REAL*8  w1(*),w2(*)
      REAL*8  xyz(3,*),xyz0(3,*),xyzfit(3,*),qt(4),error
      INTEGER ntap

      INTEGER i,iret
      REAL*8  dcm(3),q(0:7)
      CHARACTER*80 msg
      

      DO i=1,ntap
          w2(i)=w1(i)
      END DO

      CALL rigfit(0,ntap,xyz,xyz0,w1,w2,q,dcm,xyzfit,error,iret,msg)

      DO i=1,4
         qt(i)=q(i-1)
      END DO

      RETURN
      END
