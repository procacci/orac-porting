      SUBROUTINE erf_corr_cutoff(oc,delew,rkcut,nfft1,nfft2,nfft3)

************************************************************************
*   Time-stamp: <98/02/10 12:14:48 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  oc(3,3),delew,rkcut
      INTEGER nfft1,nfft2,nfft3

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  aux,mhat1,mhat2,mhat3

*----------------------- EXECUTABLE STATEMENTS ------------------------*


c--find out kcut-off in PME
c--To find out Kmax must divide by 2 nfftx, must divide
c--by 2 the OC orac matrix and must multiple by 2pi.    
      aux = 0.5*0.5*2.d0*pi
      mhat1 = aux*(oc(1,1)*nfft1+oc(2,1)*nfft2+oc(3,1)*nfft3)
      rkcut = mhat1
      mhat2 = aux*(oc(1,2)*nfft1+oc(2,2)*nfft2+oc(3,2)*nfft3)
      if(mhat2.lt.mhat1) rkcut = mhat2
      mhat3 = aux*(oc(1,3)*nfft1+oc(2,3)*nfft2+oc(3,3)*nfft3)
      if(mhat3.lt.mhat2) rkcut = mhat3
c--decrease cut-off to avoid error at the boundary 
      rkcut = 0.9*rkcut 

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
