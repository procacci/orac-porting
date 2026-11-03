      SUBROUTINE read_rem(err_args,err_unr,err_end,err_fnf)

************************************************************************
*   Time-stamp: <98/03/22 20:53:55 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Nov 18 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar
      use rem

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'mpif.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,nsevere,nwarning,j
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*15 err_fnf
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      DATA sep/' ',','/comm/'#','!'/
      REAL*8 aux
      INTEGER i
      INTEGER nproc,iproc,ierr
      real*8 dum, dum1, dum2, dum3
      CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*

      rem_groups = 0
      rem_segkind = 0
      rem_segment = .FALSE.
      rem_adjust_ratio=0.d0
      
c=======================================================================
c     Environment parser starts here 
c=======================================================================

      j=0
      nsevere = 0 
      nwarning = 0 
      line(79:80)='  '
      rem_printd=1
100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF
      
      if(mod(nproc,2).eq.0) rem_run=.true. !!!! the rem algorithm works only with an even number of replicas
      if(nproc.eq.1) rem_run=.true.        !!!! or with one replica but no exchanges are attempted

c==== Command  SCALING================================================
      
c read the extremum value of the scaling factor      
      IF(strngs(1).EQ. 'SETUP' ) THEN
        IF((nword.eq.5).or.(nword.eq.3)) THEN 
          IF(nword.eq.5) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) rem_factor_max(1)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rem_factor_max(2)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) rem_factor_max(3)
            CALL fndfmt(1,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) cdist
          ELSE                  ! if orac finds only one scaling factor, sets all the scaling factors to that value
            CALL fndfmt(2,strngs(2),fmt) 
            READ(strngs(2),fmt,err=20) rem_factor_max(1)
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) cdist
            rem_factor_max(2) = rem_factor_max(1)
            rem_factor_max(3) = rem_factor_max(1)
          ENDIF
          
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'2 or 4'
          CALL xerror(errmsg,80,1,30)
        END IF
        
c==== Command  STEP ================================================
        
c read the REM step time 
      ELSEIF(strngs(1).EQ. 'STEP' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) rem_ts
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1'
          CALL xerror(errmsg,80,1,30)
        END IF

c read the ALWAYS_ACCEPT

c==== Command  BATTERIES  ================================================
        
c read the REM step time 
      ELSEIF(strngs(1).EQ. 'BATTERIES' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) nbatteries
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  SEGMENT================================================
        
c read the first and the last atom of the segment
      ELSEIF(strngs(1).EQ. 'SEGMENT' ) THEN
        rem_segment = .TRUE.
c------- read the line
700     READ(knlist,'(a78)',END=600) line(1:78)
        CALL wrenc(kprint,line)
        IF(line(1:1) .EQ. '#') GOTO 700
        CALL parse(line,sep,2,comm,strngs,40,nword,
     x       iret,errmsg)
        IF(strngs(1) .EQ. 'define') THEN
          IF(nword.eq.3) THEN 
            rem_groups = rem_groups + 1
            if(rem_groups.gt.maxseg) THEN 
              nsevere = nsevere + 1
              errmsg=
     &  "In &REM toomany segments defined: increase maxseg in rem.f90"
              CALL xerror(errmsg,80,1,30)
            ENDIF 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) rem_group(1,rem_groups)
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rem_group(2,rem_groups)
          ELSE
            nsevere = nsevere + 1
            errmsg=err_args(3) //'2'
            CALL xerror(errmsg,80,1,30)
          END IF
        ELSEIF(strngs(1) .EQ. 'skip_force') THEN
          lskip_force=.true.  
        ELSEIF(strngs(1) .EQ. 'kind') THEN
          IF(nword.eq.2) THEN 
            IF(strngs(2) .EQ. 'intra') THEN
              rem_segkind = 1
            ELSEIF(strngs(2) .EQ. 'inter') THEN
              rem_segkind = 2
            ENDIF
! do nothing  
          ELSE
            nsevere = nsevere + 1
            errmsg=err_args(3) //'1'
            CALL xerror(errmsg,80,1,30)
          ENDIF
          
        ELSE IF(strngs(1) .EQ. 'defaults') THEN
          CONTINUE
          
        ELSE IF(strngs(1).EQ. ' ') THEN
          CONTINUE
          
        ELSE IF(strngs(1) .EQ. 'END') THEN
          GOTO 100
        ELSE
c---        could not find SUBCOMMAND of END
          errmsg=err_unr(2) // strngs(1)// ' or missing END'
          call xerror(errmsg,80,1,30)
          nsevere = nsevere + 1 
        END IF
        GOTO 700
        
        
c==== Command  PRINT================================================
        
c read the REM print time
      ELSEIF(strngs(1).EQ. 'PRINT' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) rem_print
        ELSE IF(nword.eq.3) THEN
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) rem_print
          CALL fndfmt(2,strngs(3),fmt)
          READ(strngs(3),fmt,err=20) rem_adjust_ratio
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(1) //'1 or 2'
          CALL xerror(errmsg,80,1,30)
        END IF
c read the REM_DIAGNOSTIC print time 
      ELSEIF(strngs(1).EQ. 'PRINT_DIAGNOSTIC' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) rem_printd
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  PRINT_ENERGY ===================================================
C ==

      ELSE IF(strngs(1).EQ. 'PRINT_ENERGY' ) THEN
        IF(nword == 4) then 
          prnt_remdt=.TRUE.
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) frem
          IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
              CALL openf(krem,strngs(4),'FORMATTED','OLD',0)
            ELSE
              CALL openf(krem,strngs(4),'FORMATTED','NEW',0)
            END IF
          ELSE
            errmsg=err_args(3) //'3'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
          END IF
        ELSE
          errmsg=err_args(1)//'3'
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
        ENDIF


c==== Command BLANK LINE===============================================
        
      ELSE IF(strngs(1).EQ. ' ') THEN
         GOTO 100
c==== Begininning of next ENVIRONMENT =================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF

      GO TO 100

600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_REM'
         call xerror(errmsg,80,1,2)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
         call MPI_FINALIZE(ierr)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_REM'
         call xerror(errmsg,80,1,2)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
         call MPI_FINALIZE(ierr)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_REM'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_REM'
         call xerror(errmsg,80,1,1)
      ENDIF    
      
      RETURN

c==============================================================================
c     Errors were found
c==============================================================================


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character?? '
      CALL xerror(errmsg,80,1,2)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
      call MPI_FINALIZE(ierr)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
