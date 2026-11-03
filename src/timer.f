!#ifdef Linux
      SUBROUTINE timer(vfcp,tfcp,elapse)
      REAL*8 t2(2),elapse,vfcp,tfcp,omp_get_wtime
      integer t1
      integer clock_rate, clock_max 
      tfcp=0.d0
      CALL CPU_TIME(tfcp)
#ifdef _OMP_       
      elapse=omp_get_wtime()
#else
      call system_clock ( t1, clock_rate, clock_max )
      elapse=t1/real(clock_rate)
#endif
      
      RETURN
      END
!#endif
