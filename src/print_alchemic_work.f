      subroutine print_alchemic_wrk(rtime,ladd,laddq,lrm,lrmq
     &     ,wrk_alchemy,kplot_alch0,efact)
      real*8 ladd,laddq,lrm,lrmq,wrk_alchemy,rtime,efact
      integer kplot_alch
      write(kplot_alch0,10) rtime/1000.d0,(1.d0-ladd),(1.d0-laddq),lrm
     &     ,lrmq,wrk_alchemy*efact/1000.d0
10    format(f16.3,2f10.5,3x,2f10.5,5x,f10.3) 
      return
      end
