      subroutine add_hill(nstep,time)
! add an hill to the hills list or move an old hill to a new position
! ptint out results
      
C======================= DECLARATIONS ==================================

      use unit
      use meta

      implicit none
      
C----------------------- ARGUMENTS -------------------------------------

      integer nstep
      real*8 time
      
C-------------------- LOCAL VARIABLES ----------------------------------

      integer i,ii,irc
      real*8 elapt,elapt0
      real*8 bltz,pott,meta_pot,ranf,num0,num1
      external ranf
      
C==================== EXECUTABLE STATEMENTS ============================

      call cpu_time(elapt)
      
      meta_step = meta_step + 1
      if(meta_n > metadim_nhill) then 
        write(*,*) 'actual number of hills is ', meta_n
        write(*,*) 'max number of hills is ', metadim_nhill
        write(*,*)
     &       'change metadim_nhill parameter in meta.f90 and recompile'
        stop
      end if
      pott = meta_pot(meta_rc)
      if(pott.gt.meta_maxv) meta_maxv = pott

      meta_n = meta_n + 1 
      do irc = 1,meta_nrc
        meta_list(meta_n,irc) = meta_rc(irc)
      enddo
      
      if(meta_tempered) then    ! accept a new hill with a probability exp(-V_max(t)/T')
        num0 = ranf(0.d0)
        bltz = dexp(-(meta_maxv*efact)/(gascon*meta_temperature))
        if(num0.gt.bltz) meta_n = meta_n - 1 ! if bigger, remove the last hill
      elseif(meta_wtempered) then ! accept a new hill with a probability exp(-V(s,t)/T') (well tempered algorithm)
        num0 = ranf(0.d0)
        bltz = dexp(-(pott*efact)/(gascon*meta_temperature))
        if(num0.gt.bltz) meta_n = meta_n - 1 ! if bigger, remove the last hill
      endif
      
      ii = meta_n
      
! print out results
      
      write(kmeta,'(8(f12.3,1x))') nstep*time, (meta_rc(irc), irc = 1
     &     ,meta_nrc), -meta_pot(meta_rc)*efact/1000.d0, elapt
     &     -elapt0
      
      if(metaw) then 
        if(mod(nstep*time,meta_print).lt.1.d-10) then 
          close(kmetaw)
          open(unit=kmetaw,file=metaw_file,form='FORMATTED')
          write(kmetaw,'(a1,1x,i8,6f8.3)') '#',meta_n, meta_height 
     &         *(unite * avogad) / 1000.0, (meta_width(i), i = 1
     &         ,meta_nrc)
          do i = 1,meta_n
            write(kmetaw,'(5(f12.3,1x))') (meta_list(i,irc), irc = 1
     &           ,meta_nrc)
          enddo
        endif
      endif
      
      
      elapt0 = elapt
      
      return
      end

