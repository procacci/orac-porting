      subroutine mts_plot_fragm(i1,i2,beta,xp0,yp0,zp0,ntap)
************************************************************************
*                                                                      *
*     Print snapshots of a given protein fragment                      *
*                                                                      *
************************************************************************
      use unit
      
      implicit none

      real*8  xp0(*),yp0(*),zp0(*),co(3,3),x,y,z
      CHARACTER*7 beta(*)
      integer ntap,i,i1,i2

*==================== EXECUTABLE STATEMENTS ============================
      do i = i1, i2
         x = xp0(i)
         y = yp0(i)
         z = zp0(i)
         write(kplot_fragm,1) beta(i)(1:1),x,y,z
      enddo
1     format(a1,1x,4f11.4)
*================= END OF EXECUTABLE STATEMENTS ========================
      return
      end
