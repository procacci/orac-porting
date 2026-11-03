c =====================================================================
      SUBROUTINE  steer_along_path(lstretch,nbonds_added,lbend
     &     ,nbends_added,litor,nitors_added,strbonds,strbends,strtors
     &     ,vpotbo,vpotbe,vpotit,kpotbo,kpotbe,kpotit,radfact,timerec
     &     ,ntimes,pathbo,pathbe,pathto,activebo,activebe,activeto,rtime
     &     ,eqdist,eqang,eqdied,str_vel,ben_vel,tor_vel,indexslice)
c =====================================written by P.Procacci 2009======
c  STEER_ALONG_PATH update the equilibrium steering potential and
c  calculate the actual TIMESLICE from given TIMERES for current time
c  RTIME  and calculate corresponding beginning and end of the time slice 
c  equilibrium values and velocities  
c =====================================================================

*======================== DECLARATIONS ================================*

      IMPLICIT none


*----------------------- ARGUMENTS -------------------------------------

      INTEGER lstretch,nbonds_added,lbend,nbends_added,litor
     &     ,nitors_added,ntimes,indexslice(2)
      INTEGER strbonds(*),strbends(*),strtors(*)
      REAL*8 vpotbo(*),vpotbe(*),vpotit(*),kpotbo(*),kpotbe(*),kpotit(*)
     &     ,radfact,timerec(*),rtime,eqdist(*),eqang(*),eqdied(*)
     &     ,str_vel(*),ben_vel(*),tor_vel(*),pathbo(100,*),pathbe(100,*)
     &     ,pathto(100,*)
      INTEGER activebo(100,*),activebe(100,*),activeto(100,*)  
      
*-------------------LOCAL VARIABLES ------------------------------------

      INTEGER  i,jb
      REAL*8   timeslice
      
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

c      do i=1,ntimes
c        write(6,888) rtime, timerec(i),i, indexslice(1),indexslice(2)
c     &       ,pathbo(i,1),pathbe(i,1),pathto(i,1) 
c888     format(2f10.2,3i5,3f9.3) 
c      end do
c====================================================================
c--   if time slice is found, then  we are in steering regime
c =====================================================================

*---- Bonds 
      jb=0
      if(indexslice(1).ge.1) THEN
        do i=1,nbonds_added
          if(strbonds(i).EQ.1) THEN
            jb=jb+1
*--         compute velocity in the current time slice. 
            str_vel(i) = (pathbo(indexslice(2),jb)-pathbo(indexslice(1)
     &           ,jb))/timeslice
*--         update time dependent orac potential according to current
*--         rtime 
            vpotbo(i+lstretch-nbonds_added) = pathbo(indexslice(1),jb) +
     &           str_vel(i)*(rtime-timerec(indexslice(1)))
*--         update current timeslice starting equilibrium distance 
            eqdist(i)= pathbo(indexslice(1),jb) 
c            write(6,1766) "BOND",rtime,vpotbo(i+lstretch-nbonds_added)
c     &           ,indexslice(1),jb,pathbo(indexslice(1),jb), str_vel(i),
c     &           rtime,timerec(indexslice(1))
            if(activebo(indexslice(1),jb).eq.0) THEN 
c-            driven coordinates is switched off 
              kpotbo(i+lstretch-nbonds_added) = 0.0
              str_vel(i)=0.d0 
            end if
          else
            vpotbo(i+lstretch-nbonds_added)= eqdist(i)
          end if
        end do
        
*---    Bends 
        jb=0
        do i=1,nbends_added
          if(strbends(i).eq.1) THEN
            jb=jb+1
*--         compute svelocity of the time slice. 
            ben_vel(i) = radfact*(pathbe(indexslice(2),jb)
     &           -pathbe(indexslice(1),jb))/timeslice
*--         update time dependent orac potential 
            vpotbe(i+lbend-nbends_added) = radfact*pathbe(indexslice(1)
     &           ,jb)+ben_vel(i)*(rtime-timerec(indexslice(1)))
*--         update current timeslice starting equilibrium distance 
            eqang(i)= pathbe(indexslice(1),jb) 
            if(activebe(indexslice(1),jb).eq.0) THEN 
c-            driven coordinates is switched off 
              kpotbe(i+lbend-nbends_added)= 0.0
              ben_vel(i)=0.d0
            end if
          else
            vpotbe(i+lbend-nbends_added) = eqang(i)*radfact
          end if
        end do

*---    Tors
        jb=0
        do i=1,nitors_added
          if(strtors(i).eq.1) THEN
            jb=jb+1
*--         compute svelocity of the time slice. 
            tor_vel(i) = radfact*(pathto(indexslice(2),jb)
     &           -pathto(indexslice(1),jb))/timeslice
*--         update time dependent orac potential 
            vpotit(i+litor-nitors_added) = radfact*pathto(indexslice(1)
     &           ,jb)+tor_vel(i)*(rtime-timerec(indexslice(1)))
*--         update current timeslice starting equilibrium distance 
            eqdied(i)= pathto(indexslice(1),jb) 
            if(activeto(indexslice(1),jb).eq.0) THEN 
c-            driven coordinates is switched off 
              kpotit(i+litor-nitors_added)= 0.0
              tor_vel(i)=0.d0
            end if
          else
            vpotit(i+litor-nitors_added) = eqdied(i)*radfact
          end if
        end do
      end if
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
