*======================================================================*
      subroutine steer_print(n,k,vel,lst,cord,g0,g1,t0,t1,t,w,lambda
     &     ,diff1,tstep,cost,t0r,t1r)
*======================================================================*
*  Computes work function along (harmonic) coordinate g0 for 
*  steered dynamics
*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n
      REAL*8 k(*),vel(*),g0(*),g1(*),t0,t1,t,w(*),lambda(*),t0r,t1r
     &     ,diff1(*),tstep,cord(*),cost
      INTEGER lst(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      real*8 diff0(500),xkv

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      do i=1,n
        if (lst(i).eq.1) then
          if((t.ge.t0).and.(t.le.t1)) THEN 
            diff0(i)=diff1(i)
            xkv = k(i)*vel(i)*4.184 
            lambda(i)=cost*g0(i)+vel(i)*(t-t0)
            diff1(i)= (cord(i)-lambda(i))
c--         In the funcion "work" there is a 2 factor from the 1st der. that
c--         is canceling with a 1/2 factor from trapez. numeric. integration 
            w(i) = w(i) - xkv*(diff1(i)+diff0(i))*tstep
          ELSE if(t.ge.t1r.and.t.le.t0r) THEN
            diff0(i)=diff1(i)
            xkv = k(i)*vel(i)*4.184 
            lambda(i)=cost*g1(i)-vel(i)*(t-t1r)
            diff1(i)= (cord(i)-lambda(i))
c--         In the funcion "work" there is a 2 factor from the 1st der. that
c--         is canceling with a 1/2 factor from trapez. numeric. integration 
            w(i) = w(i) - xkv*(diff1(i)+diff0(i))*tstep
          end if
        end if

      end do

      return 
      end
