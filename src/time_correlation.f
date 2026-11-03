      SUBROUTINE time_correlation(length_fft,f,g,wsave1,fcorr,fstep,w1
     &     ,w2)

************************************************************************
*   Time-stamp: <97/12/05 19:24:41 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Nov 26 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER length_fft
      REAL*8  f(*),g(*),fcorr(*),wsave1(*),fstep
      COMPLEX*16 w1(*),w2(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ns,no2,nt2,j
      REAL*8  fc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*======================================================================*
*--- Do zero padding ---------------------------------------------------
*======================================================================*

      nt2=length_fft
      no2=nt2/2
      DO i=no2+1,2*nt2
         f(i)=0.0D0
         g(i)=0.0D0
      END DO

*======================================================================*
*--- Do two transforms in one shot -------------------------------------
*======================================================================*


      CALL twofft(f,g,w1,w2,wsave1,nt2)

      DO i=1,nt2
         w2(i)=w1(i)*DCONJG(w1(i))/DFLOAT(nt2)
      END DO

*======================================================================*
*--- Finally transform back the function to t space --------------------
*======================================================================*

      CALL dcfftb(nt2,w2,wsave1)

      DO i=1,nt2
         fcorr(i)=fcorr(i)+DREAL(w2(i))
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
