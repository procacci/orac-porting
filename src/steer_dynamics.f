*======================================================================*
      subroutine steer_dynamics(m,n,str,pot,eq,eq1,vel,t,t0,t1,t0r,t1r
     &     ,ft)
*======================================================================*
*  Computes time dependent equilibrium distance/angle/torsions for steer   
*  dynamics
*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,m
      REAL*8 pot(*),vel(*),eq(*),eq1(*),t,t0,t1,t0r,t1r,ft
      integer str(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      do i=1,n
        if (str(i).eq.1) then
          if(t1r.lt.t1) THEN 
            if (t.lt.t0) then 
              pot(i+m-n) = eq(i)*ft
            else if ((t.ge.t0).and.(t.le.t1)) then
              pot(i+m-n) = eq(i)*ft+vel(i)*(t-t0)
            else if   (t.gt.t1) then
              pot(i+m-n) = eq1(i)*ft
            endif
          else
            if (t.lt.t0) then 
              pot(i+m-n) = eq(i)*ft
            else if ((t.ge.t0).and.(t.le.t1)) then
              pot(i+m-n) = eq(i)*ft+vel(i)*(t-t0)
            else if(t.gt.t1.and.t.le.t1r) then
              pot(i+m-n) = eq1(i)*ft
            else if(t.gt.t1r.and.t.lt.t0r) then 
              pot(i+m-n) = eq1(i)*ft-vel(i)*(t-t1r)
            else if(t.ge.t0r) THEN 
              pot(i+m-n)=eq(i)*ft
            end if
          endif
        endif
      enddo                

      RETURN
      END
