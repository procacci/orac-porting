      SUBROUTINE read_inout(fconf,fsave,fplot,fplot_fragm,fplot_steer,
     &     fplottemp_stn0,fplot_stn0,fplot_alch0,fplot_center,fascii
     &     ,err_open,err_args,err_end,err_unr)

************************************************************************
*   Time-stamp: <98/03/20 19:51:10 marchi>                             *
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
      use dcd, only: fdcd,kdcd,dcd_noh
#ifdef _MPI_
      use orac_mpi
#endif
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*22 err_open
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      REAL*8 fplot_center,fconf,fsave,fplot,fascii,fplot_fragm
     &     ,fplottemp_stn0,fplot_steer,fplot_stn0,fplot_alch0

#ifdef _MPI_
      INTEGER ierr,nbeg_rest_steer,nbeg_ghost
      CHARACTER*4 chari
#endif

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,j,nsevere,nwarning,stat,ilines
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      real*8   DUM
      LOGICAL  exist,plot_read_center,plot_read_fragm
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*
      kdump_out=0
      nsevere=0
      save_all_restarts=.false.
      nwarning=0
      line(79:80)='  '
      fplot_steer=0
      fplot_alch0=0
      plot_read_center=.false.
      plot_read_fragm=.false.

c=======================================================================
c     Environment parser starts here 
c=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Complex command  DUMP ===========================================

      IF(strngs(1).EQ. 'DUMP' ) THEN
         dmprnd=.TRUE.
200      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

c------- Subcommand occupy

         IF(strngs(1) .EQ. 'occupy') THEN
            occupy_space=.TRUE.

c------- Subcommand divide_records

         ELSE IF(strngs(1) .EQ. 'atom_record') THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(3)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere=nsevere+1
            END IF
            CALL fndfmt(1,strngs(2),fmt) 
            READ(strngs(2),fmt,err=20) atom_record

c------- Subcommand write

         ELSE IF(strngs(1) .EQ. 'write' .OR. strngs(1) .EQ. 'read') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fconf
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               dmpfil=strngs(4)
               INQUIRE(FILE=dmpfil,EXIST=exist)
               IF(.NOT. exist) THEN
                  errmsg=
     &                 'Auxiliary File does not exist. Abort.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 200

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 200
         
c==== Command  TRAJECTORY ===============================================

      ELSE IF(strngs(1).EQ. 'TRAJECTORY' ) THEN
        if(kprint.ne.6) THEN  
          errmsg='Invalid option for PARALLEL version; '
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
        END IF 
         dmprnd=.TRUE.
         CALL uscrpl(strngs(2),80)
         dmpfil=strngs(2)
         INQUIRE(FILE=dmpfil,EXIST=exist)
         IF(.NOT. exist) THEN
            errmsg=
     &           'Auxiliary File does not exist. Abort.'
            CALL xerror(errmsg,80,1,30)
            nsevere=nsevere+1
         END IF
         
c==== Command  RESTART  ==================================================

       ELSE IF((strngs(1).EQ. 'SAVE' .OR. strngs(1) .EQ. 'RESTART').and
     &        .nword.gt.2) THEN
         IF(strngs(2) .NE. 'CONTINUE') THEN
           CALL fndfmt(2,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) fsave
           IF(strngs(3) .EQ. 'OPEN') THEN
             CALL uscrpl(strngs(4),80)
             INQUIRE(FILE=strngs(4),EXIST=exist)
             IF(exist) THEN
               CALL openf(kdump,strngs(4),'FORMATTED','OLD',0)
             ELSE
               CALL openf(kdump,strngs(4),'FORMATTED','NEW',0)
             END IF
             restart_file=strngs(4)
             restart_out=strngs(4)
           ELSE
             errmsg=err_open
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
         ELSE
            restart_cont=.TRUE.
         END IF
         
c==== Command  RESTART  ==================================================

       ELSE IF(strngs(1) .EQ. 'RESTART'.and.nword.eq.1) THEN
300      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 300
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
         
c-------- Subcommand read
         
         IF(strngs(1) .EQ. 'read') THEN
            IF(nword .EQ. 2) THEN
               CALL uscrpl(strngs(2),80)
               INQUIRE(FILE=strngs(2),EXIST=exist)
               restart_file=strngs(2)
               IF(.NOT. exist) THEN
                  errmsg='Restart file was not found.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               ELSE
                 CALL openf(kdump,strngs(2),'FORMATTED','OLD',0)
               END IF
            ELSE
               errmsg=err_args(3)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 300

#ifdef  _MPI_ 
c-------- Subcommand read
            
          ELSE IF(strngs(1) .EQ.
     &           'read_multiple_restarts'.OR.strngs(1).EQ.'rmr') THEN
            CALL MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
            IF(nword .EQ. 3) THEN
              CALL uscrpl(strngs(2),80)
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nbeg_rest_steer
              write(chari,'(i4.4)') nbeg_rest_steer+iproc_t
cc              restart_file= '../'//TRIM(strngs(2))//chari//'.rst'
              restart_file= TRIM(strngs(2))//chari//'.rst'
              INQUIRE(FILE=restart_file,EXIST=exist)
              IF(.NOT. exist) THEN
                errmsg='Restart file'//restart_file//' was not found.'
                CALL xerror(errmsg,80,1,30)
                nsevere=nsevere+1
              ELSE
                CALL openf(kdump,restart_file,'FORMATTED','OLD',0)
              END IF
            ELSE
              errmsg=err_args(3)//'2'
              CALL xerror(errmsg,80,1,30)
              nsevere = nsevere + 1
            END IF
            GOTO 300
          ELSE IF (strngs(1).EQ.'read_multiple_pdbs'.OR.strngs(1).EQ
     &           .'rmp') THEN
            IF(nword .GE. 3) THEN
              CALL uscrpl(strngs(2),80)
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nbeg_ghost
              write(chari,'(i4.4)') nbeg_ghost+iproc_t
              pdbfile= TRIM(strngs(2))//chari//'.pdb'
              INQUIRE(FILE=pdbfile,EXIST=exist)
              IF(.NOT. exist) THEN
                errmsg='PDB ghost file'//pdbfile/
     &               /' was not found.'
                CALL xerror(errmsg,80,1,30)
                nsevere=nsevere+1
              ELSE
                CALL openf(kconf,pdbfile,'FORMATTED','OLD',0)
              END IF
            ELSE
              errmsg=err_args(3)//'3,4 or 5'
              CALL xerror(errmsg,80,1,30)
              nsevere = nsevere + 1
            END IF
            GOTO 300
#endif  
          ELSE IF(strngs(1).EQ.'read_multiple_restarts'.OR.strngs(1).EQ.
     &           'rmr') THEN
            if(kprint.eq.6) THEN  
              errmsg='Invalid option for SERIAL version; '
              CALL xerror(errmsg,80,1,30)
              nsevere=nsevere+1
            END IF 
            GOTO 300
          ELSE IF(strngs(1) .EQ. 'rmp'.OR.strngs(1).EQ
     &           .'read_multiple_pdbs') THEN
            if(kprint.eq.6) THEN  
              errmsg='Invalid option for SERIAL version; '
              CALL xerror(errmsg,80,1,30)
              nsevere=nsevere+1
            END IF 
            GOTO 300

c-------- Subcommand write
         
         ELSE IF(strngs(1) .EQ. 'write') THEN
            IF(nword .EQ. 4) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fsave
               IF(strngs(3) .EQ. 'OPEN') THEN
                  CALL uscrpl(strngs(4),80)
                  restart_out=strngs(4)
                  CALL openf(kdump_out,strngs(4),'FORMATTED','UNKNOWN'
     &                 ,0)
               ELSE IF(strngs(3) .EQ. 'SAVE_ALL_FILES') THEN
                  CALL uscrpl(strngs(4),80)
                  restart_out=strngs(4)
                  save_all_restarts=.true.
                  CALL openf(kdump_out,strngs(4),'FORMATTED','UNKNOWN'
     &                 ,0)
                ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
                END IF
            ELSE
               errmsg=err_args(3)//'3'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 300

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 300
            
         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 300

c==== Command  PLOT     ==================================================
         
       ELSE IF(strngs(1).EQ. 'PLOT' ) THEN
         IF(nword.eq.4) THEN 
           CALL fndfmt(2,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) fplot
           IF(strngs(3) .EQ. 'OPEN') THEN
             CALL uscrpl(strngs(4),80)
             INQUIRE(FILE=strngs(4),EXIST=exist)
             IF(exist) THEN
               CALL openf(kplot,strngs(4),'FORMATTED','OLD',0)
             ELSE
               CALL openf(kplot,strngs(4),'FORMATTED','NEW',0)
             END IF
           ELSE
             errmsg=err_open
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
         ELSE IF(nword .EQ. 5) THEN
           IF(strngs(2) .EQ. 'FRAGMENT') THEN
             xyz_nocell=.false.
             plot_read_fragm=.true. 
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_fragm
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_fragm,strngs(5),'FORMATTED','OLD'
     &                ,0)
               ELSE
                 CALL openf(kplot_fragm,strngs(5),'FORMATTED','NEW'
     &                ,0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'FRAGMENT_PBC') THEN
             xyz_nocell=.true.
             plot_read_fragm=.true. 
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_fragm
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_fragm,strngs(5),'FORMATTED','OLD'
     &                ,0)
               ELSE
                 CALL openf(kplot_fragm,strngs(5),'FORMATTED','NEW'
     &                ,0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'STEER') THEN
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_steer
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_steer,strngs(5),'FORMATTED','OLD'
     &                ,0)
               ELSE
                 CALL openf(kplot_steer,strngs(5),'FORMATTED','NEW'
     &                ,0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'STEER_TEMPERATURE') THEN
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplottemp_stn0
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_steer_temp,strngs(5),'FORMATTED'
     &                ,'OLD',0)
               ELSE
                 CALL openf(kplot_steer_temp,strngs(5),'FORMATTED'
     &                ,'NEW',0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'STEER_ANALYTIC') THEN
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_stn0
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_stn0,strngs(5),'FORMATTED','OLD'
     &                ,0)
               ELSE
                 CALL openf(kplot_stn0,strngs(5),'FORMATTED','NEW'
     &                ,0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'ALCHEMY') THEN
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_alch0
             if(fplot_alch0.lt.0.d0) THEN 
               print_qwork=.true. 
               fplot_alch0=dabs(fplot_alch0)
             END IF
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_alch0,strngs(5),'FORMATTED','OLD'
     &                ,0)
!                read old work    
                 ilines=0
                 stat=0
                 do while (stat.EQ.0)
                   read(kplot_alch0,*,iostat=stat) 
                   ilines=ilines+1
                 end do
                 if(ilines.le.1) THEN
                   rewind(kplot_alch0)
                   goto 102  
                 endif
                 backspace(kplot_alch0)
                 backspace(kplot_alch0)
!                read first line
                 READ(kplot_alch0,'(a78)') line(1:78)
                 IF(line(1:1) .NE. '#') THEN  
                   backspace(kplot_alch0)
                   read(kplot_alch0,*) dum,dum,dum,dum,dum,wrk_alchemy
                 end if
102              continue
               ELSE
                 CALL openf(kplot_alch0,strngs(5),'FORMATTED','NEW'
     &                ,0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE IF(strngs(2) .EQ. 'CENTER') THEN
             plot_read_center=.true. 
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) fplot_center
             IF(strngs(4) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(5),80)
               INQUIRE(FILE=strngs(5),EXIST=exist)
               IF(exist) THEN
                 CALL openf(kplot_center,strngs(5),'FORMATTED','OLD',0)
               ELSE
                 CALL openf(kplot_center,strngs(5),'FORMATTED','NEW',0)
               END IF
             ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
             END IF
           ELSE
             errmsg= err_unr(3)//strngs(1)(1:8)
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
         ELSE
           errmsg=err_args(1)//'3'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         ENDIF
         
c==== Command  EWALD_WORK =================================================

       ELSE IF(strngs(1).EQ. 'EWALD_WORK' ) THEN
         do_workq = .TRUE. 
         IF(strngs(2) .EQ. '1') THEN
           print_qwork=.true.
         ELSE IF (strngs(2) .EQ. '0') THEN
           print_qwork=.false.
         ELSE
           errmsg=err_open
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF

c==== Command  DCD     ==================================================
         
       ELSE IF(strngs(1).EQ. 'DCD' ) THEN
         IF(nword == 4 .or. nword == 5) THEN 
           CALL fndfmt(2,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) fdcd
           IF(strngs(3) .EQ. 'OPEN') THEN
             CALL uscrpl(strngs(4),80)
             INQUIRE(FILE=strngs(4),EXIST=exist)
             IF(exist) THEN
               CALL openf(kdcd,strngs(4),'UNFORMATTED','OLD',0)
             ELSE
               CALL openf(kdcd,strngs(4),'UNFORMATTED','NEW',0)
             END IF
! name for rem energies file 
             dcd_rem_file = trim(strngs(4))//'.rem'
! no hydrogens 
             if(nword == 5) then 
               if(strngs(5) == 'NOH') then 
                 dcd_noh = .true.
               else
                 errmsg= err_unr(3)//strngs(1)(1:8)
                 CALL xerror(errmsg,80,1,30)
                 nsevere = nsevere + 1
               end if
             end if
           ELSE
             errmsg=err_open
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
         ELSE
           errmsg=err_args(1)//'3'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         ENDIF
         
c==== Command  DYNAMIC =================================================

       ELSE IF(strngs(1).EQ. 'DYNAMIC' ) THEN
         DYNAMIC = .TRUE. 
         IF(strngs(2) .EQ. 'OPEN') THEN
           CALL uscrpl(strngs(3),80)
           INQUIRE(FILE=strngs(3),EXIST=exist)
           IF(exist) THEN
             CALL openf(kdynamic,strngs(3),'FORMATTED','OLD',0)
           ELSE
             CALL openf(kdynamic,strngs(3),'FORMATTED','NEW',0)
           END IF
         ELSE
           errmsg=err_open
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF

c==== Command  ASCII ===================================================
C ==

       ELSE IF(strngs(1).EQ. 'ASCII_OUTBOX' ) THEN
         ascii_nocell=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
           CALL uscrpl(strngs(4),80)
           INQUIRE(FILE=strngs(4),EXIST=exist)
           IF(exist) THEN
             CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
           ELSE
             CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
           END IF
         ELSE
           errmsg=err_open
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
         
c==== Command  ASCII ===================================================
C ==

       ELSE IF(strngs(1).EQ. 'ASCII' ) THEN
         ascii_nocell=.FALSE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
           CALL uscrpl(strngs(4),80)
           INQUIRE(FILE=strngs(4),EXIST=exist)
           IF(exist) THEN
             CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
           ELSE
             CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
           END IF
         ELSE
           errmsg=err_open
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF

c==== Command  ASCII ===================================================
C ==

       ELSE IF(strngs(1).EQ. 'PLOT_CENTER' ) THEN
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
           CALL uscrpl(strngs(4),80)
           INQUIRE(FILE=strngs(4),EXIST=exist)
           IF(exist) THEN
             CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
           ELSE
             CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
           END IF
         ELSE
           errmsg=err_open
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
         

c==== Blank Line =======================================================
C ==


       ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE

       ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1)//strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== &END keyword =====================================================
C ==

       ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

c==== No &END found or wrong keyword ===================================
C ==
       ELSE
         errmsg= err_unr(1)// strngs(1)(1:8)// err_end
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
         call int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE EXECUTING READ_INOUT'
         CALL xerror(errmsg,80,1,2)
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_INOUT'
         call xerror(errmsg,80,1,2)
         STOP
      END IF

      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//'WARNINGS WHILE EXECUTING READ_INOUT'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_INOUT'
         call xerror(errmsg,80,1,1)
      ENDIF    
      if(kdump_out.eq.0) kdump_out = kdump
      RETURN


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? Tab Character?? ' 
      CALL xerror(errmsg,80,1,2)
      STOP
      END


