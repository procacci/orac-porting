#include "config.H"	
module spme
  INTEGER :: mfft1,mfft2,mfft3,maxord,maxt,mth,max_atm,maxn
  PARAMETER(mfft1=_FFT1_+1, mfft2=_FFT2_+1, mfft3=_FFT3_+1)
  PARAMETER( &
       maxn=_MAX_GRID_+1,maxord=_MORD_, &
       maxt=2*mfft1*mfft2*mfft3,max_atm=_SIT_SOLU_, &
       mth=maxord*max_atm)
#ifdef CRAY
  REAL(8) :: bsp_mod1(mfft1),bsp_mod2(mfft2),bsp_mod3(mfft3), &
       fftable(2*(mfft1+mfft2+mfft3+50)),ffwork(2*maxt)
#endif
#ifdef FFTW               
  REAL(8) :: bsp_mod1(mfft1),bsp_mod2(mfft2),bsp_mod3(mfft3)
  INTEGER*8 planf,planr
#else
  REAL(8) :: bsp_mod1(mfft1),bsp_mod2(mfft2),bsp_mod3(mfft3), &
       fftable(3*(maxn*4+15)),ffwork(2*maxn)
#endif
  INTEGER :: sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
end module spme








