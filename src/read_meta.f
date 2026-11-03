      SUBROUTINE read_meta(err_args,err_unr,err_end)

c***********************************************************************
c                                                                      *
c              Author:  P. Procacci                                    *
c              CECAM-ENS Lyon                                          *
c                                                                      *
c              Fri Jan 17 1996                                         *
c                                                                      *
c***********************************************************************

c--- This subroutine is part of the program ORAC ----*


c======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar
      use meta

      IMPLICIT none

c----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end
      CHARACTER*27 err_unr(4)

c------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,nsevere,nwarning,j
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*15 err_fnf
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      DATA sep/' ',','/comm/'#','!'/
      REAL*8 aux

c----------------------- EXECUTABLE STATEMENTS ------------------------*

      nwarning = 0

      line(79:80)='  '
      
      meta_run  = .TRUE. ! metadynamics run
            
c=======================================================================
c     Environment parser starts here 
c=======================================================================

      j=0
      nsevere = 0 
      nwarning = 0 
      line(79:80)='  '
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


c==== Command  RATE====================================================

      IF(strngs(1).EQ. 'RATE' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_ts
          meta_analysis = .TRUE.
        ELSEIF(nword.eq.3) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_ts
          CALL fndfmt(2,strngs(3),fmt)
          READ(strngs(3),fmt,err=20) meta_height
          IF(meta_height.lt.1.d-5) meta_analysis = .TRUE.
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1 or 2'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  TEMPERED====================================================

      ELSEIF(strngs(1).EQ. 'TEMPERED' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_temperature
          meta_tempered = .TRUE.
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  WTEMPERED====================================================

      ELSEIF(strngs(1).EQ. 'WTEMPERED' ) THEN
        IF(nword.eq.2) THEN 
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_temperature
          meta_wtempered = .TRUE.
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'1'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  ADD_BOND================================================
      
      ELSEIF(strngs(1).EQ. 'ADD_BOND' ) THEN
        meta_nrc = meta_nrc + 1
        meta_rctype(meta_nrc)='bond'
        IF(nword.eq.4) THEN 
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_atoms(1,meta_nrc)
          CALL fndfmt(1,strngs(3),fmt)
          READ(strngs(3),fmt,err=20) meta_atoms(2,meta_nrc)
          CALL fndfmt(2,strngs(4),fmt)
          READ(strngs(4),fmt,err=20) meta_width(meta_nrc)
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'3'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  ADD_BEND================================================
      
      ELSEIF(strngs(1).EQ. 'ADD_BEND' ) THEN
        meta_nrc = meta_nrc + 1
        meta_rctype(meta_nrc)='bend'
        IF(nword.eq.5) THEN 
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_atoms(1,meta_nrc)
          CALL fndfmt(1,strngs(3),fmt)
          READ(strngs(3),fmt,err=20) meta_atoms(2,meta_nrc)
          CALL fndfmt(1,strngs(4),fmt)
          READ(strngs(4),fmt,err=20) meta_atoms(3,meta_nrc)
          CALL fndfmt(2,strngs(5),fmt)
          READ(strngs(5),fmt,err=20) meta_width(meta_nrc)
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'4'
          CALL xerror(errmsg,80,1,30)
        END IF


c==== Command  ADD_TORS================================================
      
      ELSEIF(strngs(1).EQ. 'ADD_TORS' ) THEN
        meta_nrc = meta_nrc + 1
        meta_rctype(meta_nrc)='tors'
        IF(nword.eq.6) THEN 
          CALL fndfmt(1,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_atoms(1,meta_nrc)
          CALL fndfmt(1,strngs(3),fmt)
          READ(strngs(3),fmt,err=20) meta_atoms(2,meta_nrc)
          CALL fndfmt(1,strngs(4),fmt)
          READ(strngs(4),fmt,err=20) meta_atoms(3,meta_nrc)
          CALL fndfmt(1,strngs(5),fmt)
          READ(strngs(5),fmt,err=20) meta_atoms(4,meta_nrc)
          CALL fndfmt(2,strngs(6),fmt)
          READ(strngs(6),fmt,err=20) meta_width(meta_nrc)
        ELSE
          nsevere = nsevere + 1
          errmsg=err_args(3) //'5'
          CALL xerror(errmsg,80,1,30)
        END IF

c==== Command  SAVE  ===============================================
        
      ELSEIF(strngs(1) .EQ. 'SAVE') THEN
        IF(nword .EQ. 3) THEN
          metaw = .TRUE.
          CALL fndfmt(2,strngs(2),fmt)
          READ(strngs(2),fmt,err=20) meta_print
          CALL uscrpl(strngs(3),80)
          metaw_file = strngs(3)
          CALL openf(kmetaw,metaw_file,'FORMATTED','UNKNOWN',0)
        ELSE
          errmsg=err_args(3) // '2'
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
        END IF
        
c==== Command  READ  ===============================================
        
      ELSEIF(strngs(1) .EQ. 'READ') THEN
        IF(nword .EQ. 2) THEN
          metar = .TRUE.
          CALL uscrpl(strngs(2),80)
          INQUIRE(FILE=strngs(2),EXIST=exist)
          IF(.NOT.exist) THEN
            errmsg='Metadynamics restart file was not found.'
            CALL xerror(errmsg,80,1,30)
            nsevere=nsevere+1
          ELSE
            metar_file = strngs(2)
            CALL openf(kmetar,metar_file,'FORMATTED','OLD',0)
          ENDIF
        ELSE
          errmsg=err_args(3) // '1'
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
        END IF
        
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
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_HISTORY'
         call xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_HISTORY'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_HISTORY'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_HISTORY'
         call xerror(errmsg,80,1,1)
      ENDIF    
      
      RETURN

c==============================================================================
c     Errors were found
c==============================================================================


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END

