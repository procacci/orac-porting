      SUBROUTINE twofft(data1,data2,fft1,fft2,wsave,n)

************************************************************************
*   Time-stamp: <97/12/15 13:52:47 marchi>                             *
*                                                                      *
*   This routine is based on its numerical recipes equivalent          *
*   The call to the complex FFT has been changed. Also the routine     *
*   has been rewritten in double precision                             *
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

      INTEGER n
      REAL*8 data1(*),data2(*),wsave(*)
      COMPLEX*16 fft1(*),fft2(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER j,n2
      COMPLEX*16 h1,h2,c1,c2

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      c1=DCMPLX(0.5,0.0)
      c2=DCMPLX(0.0,-0.5)
      DO j=1,n
        fft1(j)=DCMPLX(data1(j),data2(j))
      END DO
      CALL dcfftf(n,fft1,wsave)
      
      fft2(1)=DCMPLX(DIMAG(fft1(1)),0.0D0)
      fft1(1)=DCMPLX(DREAL(fft1(1)),0.0D0)
      n2=n+2
      DO j=2,n/2+1
        h1=c1*(fft1(j)+DCONJG(fft1(n2-j)))
        h2=c2*(fft1(j)-DCONJG(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j)=DCONJG(h1)
        fft2(j)=h2
        fft2(n2-j)=DCONJG(h2)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
