      SUBROUTINE read_simulation(fscale,err_args,err_unr,err_open
     &     ,err_end)

************************************************************************
*   Time-stamp: <98/03/05 20:33:55 marchi>                             *
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

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret,nsevere,nwarning
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*22 err_open

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  fscale  
      INTEGER nword,i,n,j
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL exist
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      nsevere = 0 
      nwarning = 0 
      line(79:80)='  '

c=======================================================================
c     Environment parser starts here 
c=======================================================================

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

c==== Command  CONST_TEMP=============================================

      IF(strngs(1).EQ. 'CONST_TEMP' ) THEN
         hoover=.TRUE.
         DO i=2,nword,2
            IF(strngs(i) .EQ. 'OMEGA') THEN
               CALL fndfmt(2,strngs(i+1),fmt)
               READ(strngs(i+1),fmt,err=20) whoov
               
            ELSE IF(strngs(i) .EQ. 'NO') THEN
               CALL fndfmt(1,strngs(i+1),fmt)
               READ(strngs(i+1),fmt,err=20) nhoov
               
            ELSE IF(strngs(i) .EQ. 'DT') THEN
               CALL fndfmt(2,strngs(i+1),fmt)
               READ(strngs(i+1),fmt,err=20) dtemph
               
            ELSE
               errmsg=err_unr(3) // strngs(i)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END DO
         IF(nhoov .NE. 1) THEN
            n=1
 15         CONTINUE
            READ(knlist,'(a78)') line(1:78)
            CALL wrenc(kprint,line)
            IF(line(1:1) .EQ. '#') GOTO 15
            CALL parse(line,sep,2,comm,strngs,40,nword,
     x           iret,errmsg)
            DO i=1,nword
               CALL fndfmt(2,strngs(i),fmt)
               READ(strngs(i),fmt,err=20) mhoov(n+i)
            END DO
            n=n+nword
            IF( n .LT. nhoov) GOTO 15
         END IF
         
c==== Command  THERMOS================================================

      ELSE IF(strngs(1).EQ. 'THERMOS' ) THEN
         thermos=.TRUE.
c------- read the line
700      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 700
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'solute') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) qmass(2)
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            
         ELSE IF(strngs(1) .EQ. 'solvent') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) qmass(3)
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'cofm') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) qmass(1)
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'temp_limit') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) dtemph
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'always_scale') THEN
            always_scale_nose=.true.

         ELSE IF(strngs(1) .EQ. 'defaults') THEN
            CONTINUE

         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1) .EQ. 'END') THEN
            GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF
         GOTO 700

c==== Command  BUSSI ============================================

      ELSE IF(strngs(1).EQ. 'BUSSI' ) THEN
         bussi=.TRUE.
         IF(nword.eq.2) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) friction
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            errmsg="friction not given: taut is 0.1 (Default value)"
            CALL xerror(errmsg,80,1,30)
            nwarning = nwarning + 1
         END IF

c==== Command  TEMPERATURE============================================

      ELSE IF(strngs(1).EQ. 'TEMPERATURE' ) THEN
         IF(nword.ge.2) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) t
            IF(nword .EQ. 3) THEN
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) dtemp
            END IF
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  ANNEALING =============================================

      ELSE IF(strngs(1).EQ. 'ANNEALING' ) THEN
         annealing=.TRUE.
         IF(nword .EQ. 2) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) annealing_fact
            annealing_fact=DSQRT(annealing_fact)
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  ISEED===================================UNSUPPORTED=====

      ELSE IF(strngs(1).EQ. 'ISEED' ) THEN
         errmsg = err_unr(4)
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         CALL fndfmt(1,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) iseed
         
c==== Command  SCALE==================================================

      ELSE IF(strngs(1).EQ. 'SCALE' ) THEN
         if(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fscale
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  ANDERSEN===============================================

      ELSE IF(strngs(1).EQ. 'ANDERSEN' ) THEN
         if(nword.eq.2) THEN 
            landersen=.TRUE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nutime
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  BERENDSEN===============================UNSUPPORTED====

      ELSE IF(strngs(1).EQ. 'BERENDSEN' ) THEN
         errmsg = err_unr(4)
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         pressure=.TRUE.
         lberendsen=.TRUE.
         IF(MOD(nword-1,2) .NE. 0) THEN
            errmsg='In READ_SIMULATION: Wrong arguments'
     &           //' after subdirective BERENDSEN. Abort.'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            DO i=2,nword,2
               IF(strngs(i) .EQ. 'TAUT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) taut
                  taut=taut*1.0D-15
               ELSE IF(strngs(i) .EQ. 'TAUP') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) taup
                  taup=taup*1.0D-15
               ELSE IF(strngs(i) .EQ. 'COMPR') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) compressibility
               ELSE IF(strngs(i) .EQ. 'PRESS-EXT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) pext
               ELSE 
                  errmsg=err_unr(3) // strngs(i)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            END DO
            IF (cpress) THEN
               errmsg='Cannot use Parrinello-Rahman'
     &              // ' and Berendsen method at the same time.' 
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            IF (hoover) THEN
               errmsg= 'Cannot use'
     &              //' Hoover and Berendsen method at the same time.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF

c==== Command  STRESS=================================================

      ELSE IF(strngs(1).EQ. 'STRESS') THEN
         cpress=.TRUE.
         pressure=.TRUE.
         IF(MOD(nword-1,2) .NE. 0) THEN
            errmsg='Wong argument(s) or invalid syntax.'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            DO i=2,nword,2
               IF(strngs(i) .EQ. 'PRESS-EXT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) pext
               ELSE IF(strngs(i) .EQ. 'BARO-MASS') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) wpr
               ELSE IF(strngs(i) .EQ. 'COMPR') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) compressibility
               ELSE IF(strngs(i) .EQ. 'TEMP_LIMIT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) dtemppr
               ELSE 
                  errmsg=err_unr(3)//strngs(i)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            END DO
            IF(strngs(4).EQ. 'iso' ) isostress=.TRUE.
         END IF

c==== Command  ISOSTRESS==============================================

      ELSE IF(strngs(1).EQ. 'ISOSTRESS' ) THEN
         cpress=.TRUE.
         isostress=.TRUE.
         pressure=.TRUE.
         IF(MOD(nword-1,2) .NE. 0) THEN
            errmsg='Wong arguments or invalid syntax.'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            DO i=2,nword,2
               IF(strngs(i) .EQ. 'PRESS-EXT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) pext
               ELSE IF(strngs(i) .EQ. 'BARO-MASS') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) wpr
               ELSE IF(strngs(i) .EQ. 'COMPR') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) compressibility
               ELSE IF(strngs(i) .EQ. 'TEMP_LIMIT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) dtemppr
               ELSE 
                  errmsg=err_unr(3)//strngs(i)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            END DO
         END IF

c==== Command  ISOSTRESSXY=============================================

      ELSE IF(strngs(1).EQ. 'ISOSTRESSXY' ) THEN
         cpress=.TRUE.
         isostressxy=.TRUE.
         pressure=.TRUE.
         IF(MOD(nword-1,2) .NE. 0) THEN
            errmsg='Wong arguments or invalid syntax.'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            DO i=2,nword,2
               IF(strngs(i) .EQ. 'PRESS-EXT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) pext
               ELSE IF(strngs(i) .EQ. 'BARO-MASS') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) wpr
               ELSE IF(strngs(i) .EQ. 'COMPR') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) compressibility
               ELSE IF(strngs(i) .EQ. 'TEMP_LIMIT') THEN
                  CALL fndfmt(2,strngs(i+1),fmt)
                  READ(strngs(i+1),fmt,err=20) dtemppr
               ELSE 
                  errmsg=err_unr(3)//strngs(i)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            END DO
         END IF

c==== Command  SCALING ==============================================


      ELSE IF(strngs(1).EQ. 'SCALING') THEN
         IF(nword .NE. 2) THEN
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         IF(strngs(2) .EQ. 'MOLECULAR') THEN
            coupl_mol=.TRUE.
            coupl_grp=.FALSE.
            coupl_atm=.FALSE.
         ELSE IF(strngs(2) .EQ. 'GROUP') THEN
            coupl_grp=.TRUE.
            coupl_mol=.FALSE.
            coupl_atm=.FALSE.
         ELSE IF(strngs(2) .EQ. 'ATOMIC') THEN
            coupl_grp=.FALSE.
            coupl_mol=.FALSE.
            coupl_atm=.TRUE.
         ELSE
            errmsg=err_unr(3)//strngs(2)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  WRITE_PRESSURE=========================================

      ELSE IF(strngs(1).EQ. 'WRITE_PRESSURE' ) THEN
         pressure=.TRUE.
         
c==== Command  MINIMIZE===============================================

      ELSE IF(strngs(1).EQ. 'MINIMIZE' ) THEN
         if(kprint.ne.6 ) THEN 
          errmsg='Invalid option for PARALLEL version; '
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
         END IF 
         mdsim=.FALSE.
         md_respa=.FALSE.
         minimize=.TRUE.
c------- read the line
1000     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1000
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'CG') THEN
            conj_grad=.TRUE.
            steepest=.FALSE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) eps_energy
         ELSE IF(strngs(1) .EQ. 'SD') THEN
            steepest=.TRUE.
            conj_grad=.FALSE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) eps_energy
         ELSE IF(strngs(1) .EQ. 'WRITE_GRADIENT') THEN
            write_grad=.TRUE.
         ELSE IF(strngs(1) .EQ. 'AGBNP') THEN
            agbnp=.true. 
         ELSE IF(strngs(1) .EQ. 'SASA') THEN
            lsasa=.true. 
         ELSE IF(strngs(1) .EQ. 'IMPLICIT_SOLVENT') THEN
            agbnp=.true. 
          ELSE IF(strngs(1).EQ. 'DOCKING' ) THEN 
            docking=.true.
            if(nword.eq.4) THEN 
              CALL fndfmt(1,strngs(2),fmt)
              READ(strngs(2),fmt,err=20) ndocking
              CALL fndfmt(2,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) dock_radius
              CALL fndfmt(2,strngs(4),fmt)
              READ(strngs(4),fmt,err=20) vsite
            ELSE if(nword.eq.5) THEN 
              CALL fndfmt(1,strngs(2),fmt)
              READ(strngs(2),fmt,err=20) ndocking
              CALL fndfmt(2,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) dock_radius
              CALL fndfmt(2,strngs(4),fmt)
              READ(strngs(4),fmt,err=20) vsite
              CALL fndfmt(1,strngs(5),fmt)
              READ(strngs(5),fmt,err=20) jligand
            ELSE
              errmsg=err_args(1)//'4'
              CALL xerror(errmsg,80,1,30)
              nsevere=nsevere+1
            END IF
         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100

         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF
         GOTO 1000
         
         
c==== Command  FREQUENCIES ================================================

      ELSE IF(strngs(1).EQ. 'FREQUENCIES' ) THEN
         mdsim=.FALSE.
         md_respa=.FALSE.
         frequencies=.TRUE.
c------- read the line
800      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 800
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'dist_max') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) hstep_freq
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            
         ELSE IF(strngs(1) .EQ. 'no_step') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) nstep_freq
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            
         ELSE IF(strngs(1) .EQ. 'print') THEN
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kfreq,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kfreq,strngs(3),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            
         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1) .EQ. 'END') THEN
            IF(kfreq .EQ. 0) THEN
               errmsg='Use the print(FREQUENCIES) subcommand to write'
     &         //     ' frequencies and eigenvectors to file'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1 
            END IF
            GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF
         GOTO 800

c==== Command  MDSIM==================================================

      ELSE IF(strngs(1).EQ. 'MDSIM' ) THEN
         mdsim=.TRUE.
         minimize=.FALSE.

c==== Blank Line =====================================================

      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

c==== &END OR WRONG KEYWORD ==========================================

      ELSE
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         call xerror(errmsg,80,1,30)
         nsevere = nsevere + 1 
      END IF

      GOTO 100
      
600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_SIMULATION'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_SIMULATION'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_SIMULATION'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_SIMULATION'
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
