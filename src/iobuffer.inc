
#include "config.H"

*------------------ PARAMETER  STATEMENT -------------------------------

      INTEGER nbuffer,buffer_time,buffer_fft
      PARAMETER (nbuffer=_BUFFER_TIME_*_BUFFER_ATOMS_,buffer_time
     &     =_BUFFER_TIME_,buffer_fft=_BUFFER_FFT_)
      REAL*4 xpb(nbuffer),ypb(nbuffer),zpb(nbuffer)
      REAL*8 xpc(buffer_time),ypc(buffer_time),zpc(buffer_time)
     &     ,vxpc(buffer_fft),vypc(buffer_fft),vzpc(buffer_fft)
     &     ,spline_x(buffer_time),spline_y(4,buffer_time)
     &     ,xpcc(buffer_time),ypcc(buffer_time),zpcc(buffer_time)
      REAL*8 wsave1(buffer_fft*4+15)
     &     ,vacf_data(0:buffer_time-1,2),vacf_spectra(0:buffer_time-1,2)
     &     ,rms_disp(0:buffer_time-1,2),tot_rms_disp(0:buffer_time-1,2)
      COMPLEX*16 w1(buffer_fft),w2(2*buffer_fft)
