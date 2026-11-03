      SUBROUTINE windows(f,nfft,nmax)

************************************************************************
*   Time-stamp: <97/12/03 18:43:06 marchi>                             *
*                                                                      *
*  Blackman windows                                                    *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Dec  3 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 f(0:*)
      INTEGER nfft,nmax

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  pi,r,w,a0,a1,a2,a3,pi2
      INTEGER nt,i
      DATA a0,a1,a2,a3/0.40217D0,0.49703D0,0.09392D0,0.00183D0/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      pi=4.0D0*DATAN(1.0D0)
      pi2=2.0D0*pi

      nt=nmax
      DO i=1,nt
          r=pi2*DFLOAT(i+nt)/DFLOAT(2*nt)
          w=a0-a1*DCOS(r)+a2*DCOS(2.0D0*r)-a3*DCOS(3.0D0*r)
          f(i)=f(i)*w
      END DO

      DO i=nt+1,nfft
          w=0.0D0
          f(i)=f(i)*w
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
