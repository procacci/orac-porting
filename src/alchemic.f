c =====================================================================
      SUBROUTINE  alchemic(timerec,ntimes,ladded,laddedq,ladd,laddq
     &     ,lremoved,lremovedq,lrm,lrmq,rtime,indexslice)
c =====================================written by P.Procacci 2013======
c  ALCHEMY updates the lamdda's (for added or removed atoms). Alchemy  
c  calculates the actual TIMESLICE from given TIMERES for current time
c  RTIME  and calculate corresponding beginning and end of the time slice 
c  equilibrium values and velocities  
c =====================================================================

*======================== DECLARATIONS ================================*

      IMPLICIT none


*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntimes,indexslice(2)
      REAL*8 timerec(*),rtime,ladded(*),lremoved(*),ladd,lrm,laddedq(*)
     &     ,lremovedq(*),laddq,lrmq
      
*-------------------LOCAL VARIABLES ------------------------------------

      INTEGER  i,jb
      REAL*8   timeslice,lambda_added_vel,lambda_removed_vel
     &     ,lambdaq_added_vel,lambdaq_removed_vel

      
      ladd=1.d0 
      lrm=0.d0 
      laddq=1.d0 
      lrmq=0.d0 
      indexslice(1)= -1
      indexslice(2)= -1
*---  find current time slice of the path according to current time
      do i=1,ntimes-1 
        if(rtime.ge.timerec(i).and.rtime.le.timerec(i+1)) then
          indexslice(1)=i
          indexslice(2)=i+1
          timeslice=timerec(i+1)-timerec(i) 
        end if
      end do

c====================================================================
c--   if time slice is found, then  we are entered into steering regime
c =====================================================================

*---- Bonds 
      jb=0
      if(indexslice(1).ge.1) THEN
        lambda_added_vel=(ladded(indexslice(2))-ladded(indexslice(1)))
     &       /timeslice
        lambda_removed_vel=(lremoved(indexslice(2))
     &       -lremoved(indexslice(1)))/timeslice
        lambdaq_added_vel=(laddedq(indexslice(2))
     &       -laddedq(indexslice(1)))/timeslice
        lambdaq_removed_vel=(lremovedq(indexslice(2))
     &       -lremovedq(indexslice(1)))/timeslice
c--     compute the current value of tghe lamda's according to the time protocol        
        ladd=ladded(indexslice(1)) + lambda_added_vel*(rtime
     &       -timerec(indexslice(1)))
        lrm=lremoved(indexslice(1)) + lambda_removed_vel*(rtime
     &       -timerec(indexslice(1)))
        laddq=laddedq(indexslice(1)) + lambdaq_added_vel*(rtime
     &       -timerec(indexslice(1)))
        lrmq=lremovedq(indexslice(1)) + lambdaq_removed_vel*(rtime
     &       -timerec(indexslice(1)))
      end if
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      if(rtime.gt.timerec(ntimes)) THEN 
        ladd=ladded(ntimes) 
        lrm=lremoved(ntimes)
        laddq=laddedq(ntimes)
        lrmq=lremovedq(ntimes) 
      END IF
      RETURN
      END
