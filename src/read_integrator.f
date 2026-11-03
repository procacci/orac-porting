      SUBROUTINE read_integrator(err_args,err_unr,err_end)

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

      IMPLICIT none

c----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*80 errmsg
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

c------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  rtol(3),rshnb(3),rnei(3)
      INTEGER nword,i,intra(2),nonbonds(3),pme_shell,ninter,nintra
     &     ,nsevere,nwarning,j
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist,dummy1,dummy2,lpme
      DATA sep/' ',','/comm/'#','!'/

c----------------------- EXECUTABLE STATEMENTS ------------------------*


      dummy1=.FALSE.
      dummy2=.FALSE.
      lpme=.false.
      ninter=0
      nintra=0
      rtol(1) = -1.0
      rtol(2) = -1.0
      rtol(3) = -1.0
      rnei(1) = -1.0
      rnei(2) = -1.0
      rnei(3) = -1.0
      nsevere = 0
      nwarning = 0

      line(79:80)='  '

c=======================================================================
c     Environment parser starts here 
c=======================================================================

 100  READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Command  USER_INTEGRATOR =========================================
 
      IF(strngs(1).EQ. 'USER_INTEGRATOR' ) THEN
         duitgl=.TRUE.
         IF(nword .GT. 1) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) dunitgi
            IF(dunitgi .GT. 10) THEN
               errmsg= err_args(2) // '10'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) dunitgr
            IF(dunitgr .GT. 10) THEN
               errmsg=err_args(2)// '10'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            IF(dunitgi .NE. 0) THEN
               DO i=1,dunitgi
                  CALL fndfmt(1,strngs(i+3),fmt)
                  READ(strngs(i+3),fmt,err=20) duitgi(i)
               END DO
            END IF
            IF(dunitgr .NE. 0) THEN
               DO i=1,dunitgr
                  CALL fndfmt(1,strngs(i+dunitgi+3),fmt)
                  READ(strngs(i+dunitgi+3),fmt,err=20) duitgr(i)
               END DO
            END IF
         END IF

c==== Command  TIMESTEP ================================================

      ELSE IF(strngs(1).EQ. 'TIMESTEP' ) THEN
         if(nword.ne.2) THEN 
            errmsg=err_args(1)//'1'
            nsevere=nsevere+1
            call xerror(errmsg,80,1,30)
         ELSE
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) time
         ENDIF

c==== Structured Command  MTS_RESPA =====================================

      ELSE IF(strngs(1).EQ. 'MTS_RESPA' ) THEN
         md_respa=.TRUE.
c------- read the line
 155     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 155
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         
c------  subcommand "dirty"---------------------------------------------- 
         IF(strngs(1).EQ.'dirty') THEN
            clean=.false. 
            
c------  subcommand "no_adjust_itors"---------------------------------------------- 
         ELSE IF(strngs(1).EQ.'no_adjust_itors') THEN
            moveitors=.false.
            
c------  subcommand "very_cold_start"------------------------------------ 
         ELSE IF(strngs(1).EQ.'very_cold_start') THEN
            start_conf=.true. 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) distmax

c------  subcommand "energy_then_die"------------------------------------ 
         ELSE IF(strngs(1).EQ.'energy_then_die') THEN 
            energy_then_die=.true. 
            
c--------subcommand "step"----------------------------------------------- 
c        SYNOPSIS "step" type n [r [hr [dr ["pme"]]]]
c        type can be "intra" or "nonbond"  
c        n : integer - respa integers  
c        r hr dr : reals - shell radius, healing and list lenghts rs 
c        if "pme" is given pme is used
c------------------------------------------------------------------------ 
         ELSE IF (strngs(1).EQ."step") THEN 
            IF(nword.ge.3) THEN 
               IF(strngs(2).eq."intra") THEN 
                  nintra = nintra + 1  
                  if(nintra.gt.2)  THEN 
                     errmsg='Too many intra shells defined - max is 2' 
                     call xerror(errmsg,80,1,30)
                     nsevere=nsevere+1
                     go to 155
                  ENDIF
                  CALL fndfmt(1,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) intra(nintra) 
               ELSE IF (strngs(2).eq."nonbond") THEN 
                  IF(nword.ge.4) THEN 
                     ninter = ninter + 1  
                     if(ninter.gt.3)  THEN 
                        errmsg
     &                  ='Too many nonbonded shells defined - max is 3'
                        call xerror(errmsg,80,1,30)
                        nsevere=nsevere+1
                        GOTO 155
                     ENDIF
                     CALL fndfmt(1,strngs(3),fmt)
                     READ(strngs(3),fmt,err=20) nonbonds(ninter) 
                     if(strngs(4).eq.'reciprocal') THEN 
                        if(.not.lpme) THEN 
                           pme_shell=ninter
                           lpme=.true.
                        else
                           errmsg='Pme shell already assigned'
                           call xerror(errmsg,80,1,30)
                           nsevere=nsevere+1
                        endif   
                     ELSE   
                        CALL fndfmt(2,strngs(4),fmt)
                        READ(strngs(4),fmt,err=20) rshnb(ninter)
                     ENDIF   
                     IF(nword.ge.5) THEN 
                        IF(strngs(5).eq.'reciprocal') THEN 
                           if (.not.lpme) THEN 
                              pme_shell=ninter
                              lpme=.true.
                           else
                              errmsg='Pme shell already assigned'
                              call xerror(errmsg,80,1,30)
                              nsevere=nsevere+1
                           end if   
                        ELSE   
                           CALL fndfmt(2,strngs(5),fmt)
                           READ(strngs(5),fmt,err=20) rtol(ninter)
                           if(nword.ge.6) THEN
                              IF(strngs(6).eq.'reciprocal') THEN 
                                 if(.not.lpme) THEN 
                                    lpme=.true.
                                    pme_shell=ninter
                                 else
                                    errmsg='Pme shell already assigned'
                                    call xerror(errmsg,80,1,30)
                                    nsevere=nsevere+1
                                 endif   
                              ELSE   
                                 CALL fndfmt(2,strngs(6),fmt)
                                 READ(strngs(6),fmt,err=20)
     &                                rnei(ninter)
                              ENDIF
                              IF(nword.eq.7) THEN 
                                 if(strngs(7).eq.'reciprocal') THEN 
                                    if(.not.lpme) then
                                       lpme=.true.
                                       pme_shell=ninter
                                    else
                                       errmsg
     &                                 ='Pme shell already assigned'
                                       call xerror(errmsg,80,1,30)
                                       nsevere=nsevere+1
                                    end if   
                                 ELSE
                                    errmsg=err_args(2) //'6'
                                    nsevere=nsevere+1
                                    call xerror(errmsg,80,1,30)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     END IF
                  ELSE
c---                 4-th argument must be the shell radius
                     errmsg=err_args(1) // '4'
                     call xerror(errmsg,80,1,30)
                     nsevere=nsevere+1
                  END IF
               ELSE
c---              2nd argument must be either "intra" or "nonbond"
                  errmsg= err_unr(3) // strngs(2) 
                  call xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               ENDIF   
            ELSE
c---           less then 3 words is a mistake: missing argument!! 
               errmsg=err_args(1)//'3'
               call xerror(errmsg,80,1,30)
               nsevere=nsevere+1
            ENDIF
c--------subcommand "k-ewald"------------------------------------------ 
         ELSE IF(strngs(1).EQ. 'k-ewald' ) THEN
            if(nword.ne.3) THEN
               nsevere = nsevere+1
               errmsg=err_args(1)//'2'
               call xerror(errmsg,80,1,30)
            else   
               clewld=.TRUE.
               grpcut=.FALSE.
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) kl
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) km
            end if
c--------subcommand "test-times"--------------------------------------- 
         ELSE IF(strngs(1).EQ. 'test_times' ) THEN
            ltest_times=.true.
            IF(strngs(2).EQ.'OPEN') THEN
               IF(nword.eq.4) THEN 
                 CALL uscrpl(strngs(3),80)
                 INQUIRE(FILE=strngs(3),EXIST=exist)
                 IF(exist) THEN
                   CALL openf(ktest,strngs(3),'FORMATTED','OLD',0)
                 ELSE
                   CALL openf(ktest,strngs(3),'FORMATTED','NEW',0)
                 END IF
                 CALL fndfmt(1,strngs(5),fmt)
                 READ(strngs(4),fmt,err=20) itprint
               ELSE IF(nword.eq.5) THEN   
                  CALL uscrpl(strngs(3),80)
                  INQUIRE(FILE=strngs(3),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(ktest,strngs(3),'FORMATTED','OLD',0)
                  ELSE
                     CALL openf(ktest,strngs(3),'FORMATTED','NEW',0)
                  END IF
                  CALL fndfmt(1,strngs(5),fmt)
                  READ(strngs(4),fmt,err=20) itprint
                  CALL fndfmt(1,strngs(5),fmt)
                  READ(strngs(5),fmt,err=20) ifprint
               ELSE IF (nword.EQ.3) THEN   ! back-compatibility
                  CALL uscrpl(strngs(3),80)
                  INQUIRE(FILE=strngs(3),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(ktest,strngs(3),'FORMATTED','OLD',0)
                  ELSE
                     CALL openf(ktest,strngs(3),'FORMATTED','NEW',0)
                  END IF
                ELSE
                 nsevere = nsevere+1
                 errmsg=err_args(1)//'2'
                 call xerror(errmsg,80,1,30)
                 errmsg=
     &          'Syntax: test_times OPEN filename istride [jstride]'
                 call xerror(errmsg,80,1,25)
              END IF
            END IF   

c--------subcommand "p_test"-----------------------------UNSUPPORTED--- 
         ELSE IF(strngs(1).EQ. 'p_test' ) THEN
            ltest_times=.true.
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ntmtsp(1)
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ntmtsp(2)
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ntmtsp(3)
            CALL fndfmt(1,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) ntmtsp(4)
            nwarning = nwarning + 1
            errmsg = err_unr(4) // strngs(1)
            call xerror(errmsg,80,1,11)
c--------subcommand "s_test"------------------------------UNSUPPORTED-- 
         ELSE IF(strngs(1).EQ. 's_test' ) THEN
            ltest_times=.true.
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ntmtss(1)
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ntmtss(2)
            CALL fndfmt(1,strngs(4),fmt)
            nwarning = nwarning + 1
            errmsg = err_unr(4) // strngs(1)
            call xerror(errmsg,80,1,11)
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
         GOTO 155
         
      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
      ELSE
c---     could not fine COMMAND of &END
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         call xerror(errmsg,80,1,30)
         nsevere = nsevere + 1 
      END IF
      GOTO 100
      
 600  CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j)//' ERRORS WHILE EXECUTING READ_INTEGRATOR'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_INTEGRATOR'
         call xerror(errmsg,80,1,2)
         STOP
      END IF

c==============================================================================
c     input verification part start here 
c==============================================================================

      if(.not.md_respa) go to 400
c---  find pme shell-----------------------------------------------------------    
      if(lpme) THEN 
         grpcut=.FALSE.
         if(pme_shell.eq.1) shell_pme='m'
         if(pme_shell.eq.2) shell_pme='l'
         if(pme_shell.eq.3) shell_pme='h'
      endif 
c---  sets INTRA subdivision n0respa and n1respa according to input-----------
      if(nintra.eq.2) THEN  
         n0respa=intra(1) 
         n1respa=intra(2) 
      else if (nintra.eq.1) THEN 
         n0respa=1
         n1respa=intra(1)
      else   
         n0respa=1
         n1respa=1
      endif

c---  sets NONBONDED subdivision m l h according to input---------------------
      
c     IF (ninter.eq.1) only one shell is defined
      IF (ninter.eq.1) THEN 
         mrespa = 1
         lrespa = 1 
         if(nonbonds(1).ne.1) THEN 
            errmsg='third argument of last "step" keyword set to 1'
            call xerror(errmsg,80,1,11)
            nwarning=nwarning + 1
         end if
         rcutm = rshnb(1)
         if (rtol(1).lt.0) THEN 
            rtolm=0.3
         ELSE
            rtolm=rtol(1)
         ENDIF   
         if (rnei(1).lt.0) THEN 
            rneim = 0.3
            rneih = rtolm+1.5
         ELSE
            rneim = 0.0
            rneih = rnei(1)+rtolm
         ENDIF   
         rcutl = rcutm 
         rtoll = 0.0
         rneil = 0.6
         rcuth = rcutm 
         rtolh = 0.0
      ENDIF
c     IF (ninter.eq.2) only two shell is defined
      IF (ninter.eq.2) THEN 
         mrespa = nonbonds(1)
         lrespa = 1 
         if(nonbonds(2).ne.1) THEN 
            errmsg='third argument of last "step" keyword set to 1'
            call xerror(errmsg,80,1,11)
            nwarning=nwarning + 1
         ENDIF
         rcutm = rshnb(1)
         if (rtol(1).lt.0) THEN 
            rtolm=0.3
         ELSE
            rtolm=rtol(1)
         ENDIF   
         if (rnei(1).lt.0) THEN 
            rneim=0.4
         ELSE
            rneim = rnei(1)
         ENDIF   
         rcutl = rshnb(2)
         rcuth = rcutl
         rtolh=0.0
         if (rtol(2).lt.0) THEN 
            rtoll=0.3
         ELSE
            rtoll=rtol(2)
         ENDIF   
         if (rnei(2).lt.0) THEN 
            rneil= 0.4
            rneih= 1.8
         ELSE
            rneil= 0.4
            rneih = rnei(2) + rtoll
         ENDIF   
         
      ENDIF
c     IF (ninter.eq.3) three shell is defined
      IF (ninter.eq.3) THEN 
         mrespa = nonbonds(1)
         lrespa = nonbonds(2) 
         if(nonbonds(3).ne.1) then 
            errmsg='third argument of last "step" keyword set to 1'
            call xerror(errmsg,80,1,11)
            nwarning=nwarning + 1
         end if
         rcutm = rshnb(1)
         if (rtol(1).lt.0) THEN 
            rtolm=0.3
         ELSE
            rtolm=rtol(1)
         ENDIF   
         if (rnei(1).lt.0) THEN 
            rneim=0.4
         ELSE
            rneim = rnei(1)
         ENDIF   
         rcutl = rshnb(2)
         if (rtol(2).lt.0) THEN 
            rtoll=0.3
         ELSE
            rtoll=rtol(2)
         ENDIF   
         if (rnei(2).lt.0) THEN 
            rneil= 0.5
         ELSE
            rneil= rnei(2)
         ENDIF   
         rcuth = rshnb(3)
         if (rtol(3).lt.0) THEN 
            rtolh=0.3
         ELSE
            rtolh=rtol(3)
         ENDIF   
         if (rnei(3).lt.0) THEN 
            rneih= 1.5
         ELSE
            rneih= rnei(3)
         ENDIF   
      ENDIF
 400  CONTINUE
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//'WARNINGS WHILE EXECUTING READ_INTEGRATOR'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_INTEGRATOR'
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

c----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
