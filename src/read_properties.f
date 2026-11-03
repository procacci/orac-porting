      SUBROUTINE read_properties(fxrms,favg,favg_xrms,fvi,gofr_fprint
     &     ,ftop_print,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf
     &     ,fnovaf,fdipole,ffragm_dist,fhbonds,frms,fgyr,err_open
     &     ,err_args,err_unr,err_end)

************************************************************************
*   Time-stamp: <98/03/20 15:37:40 marchi>                             *
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
      
      INTEGER iret
      CHARACTER*80 errmsg
      REAL*8  fxrms,fvi,gofr_fprint,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf
     &     ,fdipole,fnovaf,favg,favg_xrms,ftop_print
     &     ,ffragm_dist,fhbonds,frms,fgyr
      CHARACTER*22 err_open
      CHARACTER*37 err_args(*)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*80  itpfile

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,j,nsevere,nwarning,n_res,n_bonds
     &     ,n_bendings,n_ptors,n_itors,n_atoms,naux1,naux2
      CHARACTER*80 line,strngs(40),lined
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist,lprint
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


c=======================================================================
c     Environment parser starts here 
c=======================================================================

      nfragm = 0
      nsevere=0
      nwarning=0
      n_res=0
      n_bonds=0
      n_bendings=0
      n_ptors=0
      n_itors=0
      n_atoms=0

      DO i=1,80
         lined(i:i)=' '
      END DO
      line(79:80)='  '
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

c==== Command  TIME_CORRELATIONS =====================================

      IF(strngs(1).EQ. 'TIME_CORRELATIONS') THEN
         time_corr=.TRUE.
1200     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

c----- subcommand vacf

         IF(strngs(1) .EQ. 'vacf' ) THEN
            vacf=.TRUE.
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kvaf,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kvaf,strngs(3),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         
c---- subcommand divide_step
c---  provide number of intepolated points between data points

         ELSE IF(strngs(1).EQ. 'divide_step') THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(3)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere=nsevere+1
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) divide_spline
            END IF

c----- subcommand diffusion

         ELSE IF(strngs(1).EQ. 'diffusion') THEN
            diffusion=.TRUE.
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kdiff,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kdiff,strngs(3),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'atoms') THEN
            CALL parse_numbers(err_unr,strngs,nword,corr_atoms
     &           ,n_atoms,nsevere)

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 1200

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1200

c==== Command  PRINT_DIPOLE===========================================

      ELSE IF(strngs(1).EQ. 'PRINT_DIPOLE') THEN
         not_time_corr=.TRUE.
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fdipole
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kdipole,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kdipole,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_unr(3)//strngs(3)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command FIELD_COFACTOR==========================================

      ELSE IF(strngs(1).EQ. 'FIELD_COFACTOR') THEN 
         not_time_corr=.TRUE.
         IF(nword .LT. 2) THEN
            errmsg=' Ivalid synatx: FIELD_COFACTOR r1 [r2 r3 [OPEN fn]]'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF   
         lform_field = .false.
         IF(strngs(5) .EQ. 'OPEN') THEN
            IF(strngs(7).NE.'FORMATTED'.AND.strngs(7).NE.'UNFORMATTED')
     &          strngs(7)='UNFORMATTED' 
            IF(strngs(7).EQ.'FORMATTED') lform_field=.true. 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fvi
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) cut_field
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) alphaf
            CALL uscrpl(strngs(6),80)
            INQUIRE(FILE=strngs(6),EXIST=exist)
            IF(exist) THEN
               CALL openf(kvi,strngs(6),strngs(7),'OLD',0)
            ELSE
               CALL openf(kvi,strngs(6),strngs(7),'NEW',0)
            END IF
         ELSE IF(strngs(3) .EQ. 'OPEN') THEN 
            IF(strngs(5).NE.'FORMATTED'.AND.strngs(5).NE.'UNFORMATTED')
     &          strngs(5)='UNFORMATTED' 
            IF(strngs(5).EQ.'FORMATTED') lform_field=.true. 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fvi
            cut_field=-10.d0
            alphaf=-10.d0
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kvi,strngs(4),strngs(5),'OLD',0)
            ELSE
               CALL openf(kvi,strngs(4),strngs(5),'NEW',0)
            END IF
         ELSE   
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  DIST_FRAGMENT =====================================

      ELSE IF(strngs(1).EQ. 'DIST_FRAGMENT' ) THEN
         not_time_corr=.TRUE.
         fragm_dist=.TRUE.
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ffragm_dist
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kfragm_dist,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kfragm_dist,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_unr(3)//strngs(3)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  DEF_FRAGMMENT =========================================

      ELSE IF(strngs(1).EQ. 'DEF_FRAGMENT' ) THEN
         not_time_corr=.TRUE.
         nfragm = nfragm+1
         IF(nfragm.le.20) THEN
            IF(nword .NE. 3) THEN
               errmsg=err_args(1)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fragm(1,nfragm) 
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) fragm(2,nfragm) 
            END IF
         ELSE
            errmsg='Toomany fragments defined; Max is 20'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  GOFR===================================================

      ELSE IF(strngs(1) .EQ. 'GOFR') THEN
         not_time_corr=.TRUE.
         gofr=.TRUE.
 700     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 700
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_fprint
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kgofr_sk,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kgofr_sk,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'delta') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) delrg

         ELSE IF(strngs(1) .EQ. 'average') THEN
            gofr_avg=.TRUE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_favg

         ELSE IF(strngs(1) .EQ. 'compute') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_fcomp

         ELSE IF(strngs(1) .EQ. 'use_neighbor') THEN
            gofr_neighbor=.TRUE.

         ELSE IF(strngs(1) .EQ. 'cutoff') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_cut

         ELSE IF(strngs(1) .EQ. 'intra') THEN
            gofr_intra=.TRUE.

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 700

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 700


c==== Command  HBONDS ================================================

      ELSE IF(strngs(1).EQ. 'HBONDS' ) THEN
         not_time_corr=.TRUE.
1000     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1000
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1) .EQ. 'total') THEN
            hbonds_tot=.TRUE.

         ELSE IF(strngs(1) .EQ. 'residues') THEN
            hbonds_res=.TRUE.

         ELSE IF(strngs(1) .EQ. 'radial_cutoff') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) rcut_hb
            ELSE
               errmsg=err_args(2)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'angular_cutoff') THEN
            IF(nword .EQ. 3) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) acut_hb
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) a2cut_hb
            ELSE IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) acut_hb
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fhbonds
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(khbonds,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(khbonds,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1000

c==== Command COMPARE =================================================

         ELSE IF(strngs(1).EQ. 'STRUCTURES' ) THEN
            not_time_corr=.TRUE.
1100        READ(knlist,'(a78)',END=600) line(1:78)
            CALL wrenc(kprint,line)
            IF(line(1:1) .EQ. '#') GOTO 1100

            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
            
c==== Subcommand print structure ====================================

            IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ. 'averaged')
     &           THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_str=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) favg
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kavg,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kavg,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            
c==== Subcommand print avg_xrms =====================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ.
     &              'avg_xrms') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_str=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) favg_xrms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kavg_xrms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kavg_xrms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand print inst_xrms ====================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ.
     &              'inst_xrms') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) fxrms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kxrms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kxrms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
                  CALL add_str(strngs(5),80,'_atm',4,lined)
                  INQUIRE(FILE=lined,EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kxrms_atm,lined,'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kxrms_atm,lined,'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand ins_xrms ============================================

            ELSE IF(strngs(1).EQ. 'inst_xrms') THEN
               anxrms=.TRUE.
               DO i=2,nword
                  IF(strngs(i).EQ. 'ca') THEN
                     anxca=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'backbone') THEN
                     anxbc=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'heavy') THEN
                     anxhe=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'allatoms')
     &                    THEN
                     anxal=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'all') THEN
                     anxca=.TRUE.
                     anxbc=.TRUE.
                     anxhe=.TRUE.
                     anxal=.TRUE.
                  ELSE
                     errmsg=err_unr(3) // strngs(i)
                     CALL xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                  END IF
               END DO


c==== Subcommand select_cell=========================================

            ELSE IF(strngs(1).EQ. 'select_cell') THEN
               anxrms_cell=.TRUE.

c==== Subcommand averaged ===========================================

            ELSE IF(strngs(1).EQ. 'averaged') THEN
               IF(nword .LT. 2) THEN
                  errmsg=err_args(1)//'1'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(2) .EQ. 'ca') THEN
                  avg_ca=.TRUE.
               ELSE IF(strngs(2) .EQ. 'heavy') THEN
                  avg_he=.TRUE.
               ELSE
                  errmsg=err_unr(3)//strngs(2)//err_end(1:14)
     &                 //err_end(16:20)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand print rms ==========================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ. 'rms')
     &              THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_rms=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) frms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(krms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(krms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF



            ELSE IF(strngs(1) .EQ. ' ') THEN
               GOTO 1100
               
            ELSE IF(strngs(1).EQ. 'END' ) THEN
               GOTO 100

            ELSE
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/
     &              /err_end(16:20)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 1100

c==== Command TOPOLOGY================================================

      ELSE IF(strngs(1).EQ. 'TOPOLOGY' ) THEN
         not_time_corr=.TRUE.
         prttopl=.TRUE.
900      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 900
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ftop_print
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(ktopol,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(ktopol,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'bonds') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_bonds,n_bonds
     &           ,nsevere)
            IF(n_bonds .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'bendings') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_bendings
     &           ,n_bendings,nsevere)
            IF(n_bendings .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'P-torsions') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_ptors,n_ptors
     &           ,nsevere)
            IF(n_ptors .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'I-torsions') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_itors,n_itors
     &           ,nsevere)
            IF(n_itors .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 900

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(1)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 900
         
c==== Command FORCE_FIELD ============================================

      ELSE IF(strngs(1) .EQ. 'FORCE_FIELD') THEN
         write_ff_pars=.TRUE.
         
      ELSE IF(strngs(1) .EQ. 'WRITE_ITP') THEN
         gen_itp=.TRUE.
         itpfile='file.itp'
         CALL openf(kitp,itpfile,'FORMATTED','UNKNOWN',0)
         
c==== Command GYRATION ===============================================

      ELSE IF(strngs(1).EQ. 'GYRATION') THEN
         not_time_corr=.TRUE.
         wrtgyr=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fgyr
         IF(strngs(3).EQ. 'OPEN') THEN
            
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kgyr,strngs(4),
     &              'FORMATTED','OLD',0)
            ELSE
               CALL openf(kgyr,strngs(4),
     &              'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
c==== Begininning of next ENVIRONMENT =================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

         CONTINUE
         
c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
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
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_PROPERTIES'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_PROPERTIES'
         call xerror(errmsg,80,1,1)
      ENDIF    
      IF(iret.eq.1) THEN
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF
      top_bonds(1)=n_bonds
      top_bendings(1)=n_bendings
      top_ptors(1)=n_ptors
      top_itors(1)=n_itors
      corr_atoms(1)=n_atoms

      RETURN

c==============================================================================
c     Errors were found
c==============================================================================


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      RETURN

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
