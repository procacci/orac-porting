#define PRESSURE
      SUBROUTINE verify_input(fmaxstp,fscale,fprint,fsave,fplot,fascii
     &     ,fprop,fplot_fragm,fplot_steer,fplottemp_stn0,fplot_stn0
     &     ,fplot_alch0,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms
     &     ,fvi,favg,favg_xrms,ftop_print,gofr_fprint,gofr_favg
     &     ,gofr_fcomp,fprtvaf,ftotvaf,fnovaf,fdipole,ffragm_dist
     &     ,fhbonds,frms,fgyr,fabmd,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/07 11:44:20 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Nov 19 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar
      use spme
      use meta
      use rem
      use sge, only: sge_tsave_work, sge_tsave_pmf, sim_tempering,
     &     sim_tempering_sge, sge_ts, sge_f, sge_print_wham, sge_print
      use dcd, only: ndcd,fdcd
#ifdef  _MPI_
      use orac_mpi, only:nfile,iproc,nproc
#endif
      
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  fmaxstp,fscale,fprint,fsave,fplot,fascii,fprop,fplot_fragm
     &     ,fplot_steer,fplottemp_stn0,fplot_stn0,fplot_alch0
     &     ,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms,favg
     &     ,favg_xrms,ftop_print,gofr_fprint,gofr_favg
     &     ,gofr_fcomp,fprtvaf,ftotvaf,fnovaf,fdipole,fvi,ffragm_dist
     &     ,fhbonds,frms,fgyr,fabmd
      INTEGER iret,nwarning,nsevere
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 n_timestep,aux,cut,tsign,time0
      LOGICAL abort1,abort1a,abort1b,abort1c,abort2,near0,ok
      INTEGER i,idum1,idum2
      CHARACTER*11 integers 
      CHARACTER*2 nerr
      logical exist

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      aux=DABS(time)
      tsign=aux/time
      time=aux
      write(kprint,'(5x,a)')
     &     'Checking Input ......              ---->'
      integers=" 1234567890"
      nerr="  "
      nwarning = 0
      nsevere  = 0
      IF(.NOT. replicate) THEN
         aaxis=aaxis/DFLOAT(icl)
         baxis=baxis/DFLOAT(icm)
         caxis=caxis/DFLOAT(icn)
      END IF
         
      IF(slt_exist .AND. slv_create) THEN
         slt_create=.TRUE.
      END IF

      IF(slt_exist .AND. sgroup) THEN
         slt_create=.TRUE.
      END IF

      IF(debug.and.tpgfil) THEN 
         errmsg=
     &'DEBUG(&RUN) ignored when READ_PFR_BIN(&PARAMETERS) is'
     & / / ' specified'       
         call xerror(errmsg,80,1,21)
         nwarning   = nwarning  + 1 
      END IF  
      IF(grpcut .AND. pme) THEN
         errmsg=
     &'PME and GROUP_CUTOFF(&POTENTIAL) are mutually exclusive'
         call xerror(errmsg,80,1,20)
         nsevere  = nsevere  + 1 
      END IF

      IF(grpcut .AND. clewld) THEN
         errmsg
     &='EWALD and GROUP_CUTOFF(&POTENTIAL) are mutually exclusive'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(solven.and.slvpdb.and.linser) THEN
         errmsg=
     &'READ_PDB(&SOLVEN) AND INSERT(&SETUP) have been both specified'
         call xerror(errmsg,80,1,21)
         nwarning = nwarning + 1
      END IF
      
      IF(pme .AND. md_respa) THEN
         IF(shell_pme .EQ. '0') THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) no shell defined for PME'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         IF(nfft1 .EQ. 0) THEN
            errmsg=
     & 'In EWALD(&POTENTIAL) no grid speicified for pme.'  
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

      IF(clewld .AND. md_respa .AND. (.NOT. pme) ) THEN
         IF(shell_pme .EQ. '0') THEN
            errmsg=
     & 'In EWALD(&POTENTIAL) no shell defined'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      
      IF(pme) THEN
         IF(nfft1+1 .GT. mfft1 .OR. nfft2+1 .GT. mfft2 .OR. nfft3+1 .GT.
     &        mfft3) THEN
            errmsg=
     &'pme grid nfft1-2-3 in EWALD(&POTENTIAL) exceed physical'
     &/ /' dimensions'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

#ifndef PRESSURE
      IF(cpress .OR. pressure .OR. isostress) THEN
         errmsg=
     & 'In &SIMULATION(STRESS) should recompile with pressure option.'
     & / /' Change pressure.H'
            call xerror(errmsg,80,1,20)
         nsevere=nsevere+1
      END IF
#endif

*=======================================================================
*--- Validate steps                                                  ---
*=======================================================================

      IF(minimize) THEN
         time=1.0D0
         nrespa=1
         n0respa=1
         n1respa=1
         mrespa=1
         lrespa=1
      END IF

      IF( .NOT. analys) THEN
         maxstp=IDINT(fmaxstp/time)
         nprop=IDINT(fprop/time)
         nrject=IDINT(frject/time)
         nsave=IDINT(fsave/time)
         nupdte=IDNINT(fupdte/time)
         n_timestep=time/DFLOAT(mrespa*lrespa)
         time0 = n_timestep/float(n1respa*n0respa)
         maxrun=IDINT(fmaxrun/n_timestep)
         nprint=IDNINT(fprint/n_timestep)
         nconf=IDNINT(fconf/n_timestep)
         nplot_fragm=NINT(fplot_fragm/n_timestep)
         nplot_steer=NINT(fplot_steer/time)
         if(time0.gt.fplot_stn0) fplot_stn0=time0
         if(n_timestep.gt.fplot_alch0.and.fplot_alch0.gt.0.d0)
     &        fplot_alch0=n_timestep
         if(time0.gt.fplottemp_stn0) fplottemp_stn0=time0
         nplotstn0=NINT(fplot_stn0/time0)
         nplotalch0=NINT(fplot_alch0/n_timestep)
         nplottempstn0=NINT(fplottemp_stn0/time0)
         nplot_center=NINT(fplot_center/n_timestep)
         ndcd = idnint(fdcd/n_timestep)
         nascii=IDNINT(fascii/n_timestep)
         nplot=IDNINT(fplot/n_timestep)
         nxrms=IDNINT(fxrms/n_timestep)
         nvi=IDNINT(fvi/n_timestep)
         ndipole=IDNINT(fdipole/n_timestep)
         nfragm_dist=IDNINT(ffragm_dist/n_timestep)
         nhbonds=IDNINT(fhbonds/n_timestep)
         nrms=IDNINT(frms/n_timestep)
         ngyr=IDNINT(fgyr/n_timestep)
         navg_str=IDNINT(favg/n_timestep)
         navg_str_xrms=IDNINT(favg_xrms/n_timestep)
         ntop_print=IDNINT(ftop_print/n_timestep)
         gofr_nprint=NINT(gofr_fprint/n_timestep)
         gofr_navg=NINT(gofr_favg/n_timestep)
         gofr_ncomp=IDINT(gofr_fcomp/n_timestep)
         nabmd=IDINT(fabmd/n_timestep)
         IF(NINT(fscale/n_timestep) .NE. 0) scale=NINT(fscale/n_timestep
     &     )
         nrem = idnint(frem/n_timestep)
      ELSE
         nplot_fragm=IDNINT(fplot_fragm)
         nplot_center=IDNINT(fplot_center)
         nascii=IDNINT(fascii)
         nplot=IDNINT(fplot)
         nxrms=IDNINT(fxrms)
         nplot_fragm=IDNINT(fplot_fragm)
         nfragm_dist=IDNINT(ffragm_dist)
         ndipole=IDNINT(fdipole)
         navg_str=IDNINT(favg)
         navg_str_xrms=IDNINT(favg_xrms)
         ntop_print=IDNINT(ftop_print)
         gofr_nprint=IDNINT(gofr_fprint)
         gofr_navg=IDNINT(gofr_favg)
         gofr_ncomp=IDNINT(gofr_fcomp)
         nhbonds=IDNINT(fhbonds)
         nrms=IDNINT(frms)
         ngyr=IDNINT(fgyr)
         IF(vacf) THEN
            prtvaf=IDINT(fprtvaf/n_timestep)
            prtvaf=prtvaf/novaf
         END IF
      END IF

      if ( aaxis.EQ.500.0 .and.   
     &     baxis.EQ.500.0 .and.   
     &     caxis.EQ.500.0 .and.   
     &     alf.EQ.90.0 .and.   
     &     bet.EQ.90.0 .and.   
     &     gam.EQ.90.0.and. (.not.docking)) THEN 
        errmsg='BOX not given in &READ_SETUP'
        call xerror(errmsg,80,1,21)
        errmsg='Assumed box size is 500x500x500 Angs^3'
        call xerror(errmsg,80,1,26)
        nwarning = nwarning + 1
      end if
        
      if(docking) THEN
        aaxis=1000.d0
        baxis=1000.d0
        caxis=1000.d0
        write(kprint,1065)
1065    format("     Docking is ON: large box is assumed  ---->")  
      endif

        
      if( sim_tempering ) then

        if( steer_opt ) then
          errmsg='A SGE simulation is being performed [see SETUP(&SGE)]'
          call xerror(errmsg,80,1,20)
          errmsg = 'The command STEER(&RUN) is incompatible.'
          call xerror(errmsg,80,1,26)
          errmsg = 'The command STEER(&RUN) should be removed.'
          call xerror(errmsg,80,1,25)
          nsevere = nsevere + 1
        endif

        if( path_steer ) then
          errmsg='A SGE simulation is being performed [see SETUP(&SGE)]'
          call xerror(errmsg,80,1,20)
          errmsg = 'The command STEER_PATH(&RUN) is incompatible.'
          call xerror(errmsg,80,1,26)
          errmsg='The command STEER_PATH(&POTENTIAL) should be removed.'
          call xerror(errmsg,80,1,25)
          nsevere = nsevere + 1
        endif

        if( sim_tempering_sge .and. rem_segment ) then
          errmsg = 'A SGE simulation on a reaction coordinate'//
     &              ' is being performed [see SETUP(&SGE)].'
          call xerror(errmsg,80,1,20)
          errmsg = 'The command SEGMENT(&SGE) is incompatible.'
          call xerror(errmsg,80,1,26)
          errmsg = 'The command SEGMENT(&SGE) should be removed.'
          call xerror(errmsg,80,1,25)
          nsevere = nsevere + 1
        endif

        if( sge_ts /= time ) then
c--     sets sge_ts for ensemble transitions to an EXACT multiple of the MTS-timestep
          sge_ts = time * nint(sge_ts/time)
          if( sge_ts < time ) then
            sge_ts = time
            errmsg =
     & 'sge_ts(&SGE) for ensemble transition is less than the timestep.'
            call xerror(errmsg,80,1,21)
            errmsg = 'sge_ts set to TIMESTEP'
            call xerror(errmsg,80,1,26)
            nwarning = nwarning + 1
          endif
        endif


        if ( sge_print_wham .and. sge_f /= time ) then
c--     sets sge_f for printing SGE_WHAM file to an EXACT multiple of the MTS-timestep
          sge_f = time * nint(sge_f/time)
          if ( sge_f < time ) then
            sge_f = time
            errmsg =
     & 'sge_f(&SGE) to print SGE_WHAM file is less than the timestep.'
            call xerror(errmsg,80,1,21)
            errmsg = 'sge_f set to TIMESTEP'
            call xerror(errmsg,80,1,26)
            nwarning = nwarning + 1
          endif
        endif

        if( sge_tsave_work /= time ) then
c--     sets sge_tsave_work to an EXACT multiple of the MTS timestep 
          sge_tsave_work = time * nint(sge_tsave_work/time)
          if( sge_tsave_work < time ) then
            sge_tsave_work = time
            errmsg=
     &        'sge_tsave_work(&SGE) for SGE is less than the timestep.'
            call xerror(errmsg,80,1,21)
            errmsg = 'sge_tsave_work set to TIMESTEP'
            call xerror(errmsg,80,1,26)
            nwarning = nwarning + 1
          endif
        endif

        if( sge_tsave_pmf /= time ) then
c--     sets sge_tsave_pmf to an EXACT multiple of the MTS timestep 
          sge_tsave_pmf = time * nint(sge_tsave_pmf/time)
          IF( sge_tsave_pmf < time ) then
            sge_tsave_pmf = time * 1000
            errmsg=
     &        'sge_tsave_pmf(&SGE) for SGE is less than the timestep.'
            call xerror(errmsg,80,1,21)
            errmsg = 'sge_tsave_pmf set to 1000 x TIMESTEP'
            call xerror(errmsg,80,1,26)
            nwarning = nwarning + 1
          endif
        endif

        if( sge_print /= time ) then
c--     sets sge_print for printing to an EXACT multiple of the MTS-timestep
          sge_print = time * nint(sge_print/time)
          if( sge_print < time ) then
            sge_print = sge_tsave_pmf
            errmsg =
     & 'sge_print(&SGE) is less than the timestep.'
            call xerror(errmsg,80,1,21)
            errmsg = 'sge_print set to sge_tsave_pmf'
            call xerror(errmsg,80,1,26)
            nwarning = nwarning + 1
          endif
        endif

        if( sge_tsave_pmf < sge_tsave_work ) then
          errmsg =
     &         'Error in STEP(&SGE). The time interval used to update'
          call xerror(errmsg,80,1,20)
          errmsg =
     &         'the free energies is smaller than the time interval'
          call xerror(errmsg,80,1,26)
          errmsg = 'used to update the generalized work.'
          call xerror(errmsg,80,1,26)
          nsevere = nsevere + 1
        endif

#ifdef _MPI_
      ELSEIF(rem_run) THEN

        IF(rem_ts.NE.time) THEN
c--     sets rem_step for exchange to an EXACT multiple of the MTS
c--     -timestep 
          rem_ts = time * nint(rem_ts/time)
          IF(rem_ts.LT.time) THEN
            rem_ts=time
            errmsg=
     &         'rem_ts(&REM) for replica exchange is less than the times
     &tep ;'//'rem_ts set to TIMESTEP'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
          ENDIF
        ENDIF

        IF(rem_print.NE.time) THEN
c--     sets rem_print for printing to an EXACT multiple of the MTS
c--     -timestep 
          rem_print = time * nint(rem_print/time)
          IF(rem_print.LT.time) THEN
            rem_print=time
            errmsg=
     &         'rem_print(&REM) for printing is less than the times
     &tep ;'//'rem_print set to TIMESTEP'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
          ENDIF
        ENDIF
        if(nbatteries.gt.1) THEN
          if(mod(nproc,nbatteries).ne.0) THEN 
            if(iproc.eq.0) write(*,12754) nbatteries,nproc
12754       FORMAT( '************************************************' /
     &              '*          FATAL ERROR!!!                      *' / 
     &              '* -nbatteries is not a multiple of nproc-      *' / 
     &              '* nbatteries =',i4,' nproc= ',i4,'             *' / 
     &              '* Action: change nbatteries in &READ_REM       *' / 
     &       '************************************************' /)
            CALL MPI_Finalize(iret)
            STOP
          end if
        end if

#endif
      ENDIF

! when plotting a dcd trajectory during a rem simulation
! automatically print out rem energies with the same frequency
      
      if(rem_run .and. ndcd /= 0) then 
        INQUIRE(FILE=trim(dcd_rem_file),EXIST=exist)
        IF(exist) THEN
          CALL openf(krem_dcd,dcd_rem_file,'FORMATTED','OLD',0)
        ELSE
          CALL openf(krem_dcd,dcd_rem_file,'FORMATTED','NEW',0)
        END IF
      end if

      IF(str_tim1r.gt.0.d0.and.str_tim1r.lt.str_time1) THEN 
        errmsg=
     &'str_tim1r(&RUN) is less than str_tim1 in STEER dynamics;'/
     &   /' it must be str_tim1r > str_tim1'
        nprint = 1
        call xerror(errmsg,80,1,20)
        nsevere=nsevere+1
      ENDIF
      IF(.NOT. stoprun) THEN
         IF( nprint .EQ. 0) THEN
            errmsg=
     &'PRINT(&RUN) is 0 or less than TIMESTEP(&VARIABLES);'/ /
     &           ' NPRINT set to 1.'
            nprint = 1
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
         END IF
         IF( nprop .EQ. 0) THEN
            errmsg=
     &'PROPERTY(&RUN) is O or less than TIMESTEP(&VARIABLES)'/ /
     &'; NPROP set to 1000'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
            nprop=1000
         END IF
         IF( nupdte .EQ. 0 .AND. (.NOT. analys)) THEN
            errmsg=
     &'UPDATE(&POTENTIAL) is 0 or less than TIMESTEP'/ /
     &'(&VARIABLES); NUPDTE set to 4'
            call xerror(errmsg,80,1,21)
            nupdte=4
            nwarning = nwarning + 1 
         END IF
            
         IF(gofr) THEN
            IF(gofr_ncomp .EQ. 0) THEN
               errmsg=
     &'compute(GOFR&PROPERTIES) is 0 or less'/ / 
     &'than TIMESTEP(&VARIABLES); reset to 1'
               call xerror(errmsg,80,1,21)
               nwarning = nwarning + 1 
               gofr_ncomp = 1
            END IF
            IF(gofr_nprint .EQ. 0) THEN
               errmsg=
     &'print(GOFR&PROPERTIES) is 0 or less than '/ /
     &'TIMESTEP(&VARIABLES); reset to 1'
               call xerror(errmsg,80,1,21)
               nwarning = nwarning + 1
               gofr_nprint = 1
            END IF
         END IF
      END IF
      IF(linked_cell) THEN
         IF(nupdte_index .EQ. 0) THEN
            errmsg=' Frequency of calls to index cell list is zero. '
     &           / /'Check LINKED_CELL(&POTENTIAL). '
            CALL xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF((.NOT. bending) .AND. stretch) THEN
          errmsg='Cannot use stretching and bending at the'
     x  / /' same time. Abort.'
          CALL xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(rspoff .EQ. 0.0D0 .AND. (.NOT. md_respa) .AND. (.NOT. analys))
     &     THEN
         errmsg=
     &'Should specify a cutoff in &POTENTIAL (commands CUTOFF'
     &/ /' or GROUP_CUTOFF)'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(gofr .AND. gofr_neighbor) THEN
         IF(.NOT. md_respa) THEN
            cut=rspoff
            IF(rspcut+cut .LT. gofr_cut) THEN
               gofr_neighbor=.FALSE.
            END IF
            IF(gofr_cut .EQ. 0.0D0) gofr_cut=cut
         ELSE
            cut=rcuth+rtolh
            IF(rneih+cut .LT. gofr_cut) THEN
               gofr_neighbor=.FALSE.
            END IF
            IF(gofr_cut .EQ. 0.0D0) gofr_cut=rcuth+rtolh
         END IF
       END IF
       
       IF(gofr) THEN
         IF(NINT(gofr_cut/delrg) .GT. maxint) THEN
           errmsg=
     &'**ERROR: Dimensions of pair correlation are insufficient.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      
      IF(meta_run) THEN
        IF(meta_nrc .EQ. 0) THEN
          errmsg='**ERROR: No reaction coordinates 
     &have been defined for this metadynamics run.'
          call xerror(errmsg,80,1,20)
          nsevere = nsevere + 1 
        END IF
        IF(meta_ts.ne.time) THEN 
c--     sets deposition frequency for metadynamics to an EXACT multiple
c--     of the MTS
c--     -timestep 
          meta_ts=time*nint(meta_ts/time) 
          IF(meta_ts.lt.time) THEN
            meta_ts=time
            errmsg=
     &           'meta_ts(&META) for metadynamics is less than the times
     &tep ;'//'meta_ts set to TIMESTEP'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
          END IF
        ENDIF
        IF(meta_tempered.or.meta_wtempered) THEN
          IF(meta_temperature.LT.1.d-5) THEN
            meta_tempered = .FALSE.
            meta_wtempered = .FALSE.
            errmsg
     &           ='temperature was zero, performing a standard meta run'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
          ENDIF
        END IF
      END IF

#ifdef _MPI_
      if( rem_run ) then
        if ( cdist .eq. 0 .and. nflag(1) .eq. 0 ) then
          errmsg =
     &         'cannot restart a REM run' //
     &         ' without reading a restart file.'
          call xerror(errmsg,80,1,20)
          errmsg =
     &         'A conflict occurs between the restart option' //
     &         ' of command SETUP(&REM) and the'
          call xerror(errmsg,80,1,26)
          errmsg =
     &         'restart option of command CONTROL(&RUN)'
          call xerror(errmsg,80,1,26)
          nsevere = nsevere + 1
        elseif ( cdist .GT. 0 .AND. nflag(1) .GT. 0 ) THEN 
          errmsg =
     &         'can read a restart file' // 
     &         'only when restarting a previous REM run' 
          call xerror(errmsg,80,1,20)
          errmsg =
     &         'A conflict occurs between the restart option' //
     &         ' of command SETUP(&REM) and the'
          call xerror(errmsg,80,1,26)
          errmsg =
     &         'restart option of command CONTROL(&RUN)'
          call xerror(errmsg,80,1,26)
          nsevere = nsevere + 1
        endif
      endif
#endif

      if( sim_tempering ) then
        if ( cdist .eq. 0 .and. nflag(1) .eq. 0 ) then
          errmsg =
     &         'Cannot restart a SGE run' //
     &         ' without reading a restart file.'
          call xerror(errmsg,80,1,20)
          errmsg =
     &         'A conflict occurs between the restart option' //
     &         ' of command SETUP(&SGE) and the'
          call xerror(errmsg,80,1,26)
          errmsg =
     &         'restart option of command CONTROL(&RUN)'
          call xerror(errmsg,80,1,26)
          nsevere = nsevere + 1
        elseif( cdist .gt. 0 .and. nflag(1) .gt. 0) then 
          errmsg =
     &         'Can read a restart file' //
     &         ' only when restarting a previous SGE run.'
          call xerror(errmsg,80,1,20)
          errmsg =
     &         'A conflict occurs between the restart option' //
     &         ' of command SETUP(&SGE) and the'
          call xerror(errmsg,80,1,26)
          errmsg =
     &         'restart option of command CONTROL(&RUN)'
          call xerror(errmsg,80,1,26)
          nsevere = nsevere + 1
        endif
      endif


*=======================================================================
*----- Consistency checks ---------------------------------------------
*=======================================================================

      IF(minimize) THEN
        IF(near0(rspoff)) THEN
          errmsg=
     &         'Use CUTOFF(&POTENTIAL) to assign a direct space'/ /
     &         ' cutoff.'
          call xerror(errmsg,80,1,20)
          nsevere = nsevere + 1 
        END IF
        rcuth=rspoff-rspon
        rtolh=rspon
        rtoll=rcuth
        rcutl=rspoff-rspon
        rtoll=0.0
        rneih=rspcut+rspoff
        rneil=0.0
        IF(md_respa) THEN
          errmsg=
     &         'Warning MTS_STEP(&INTEGRATOR) option is'/ /
     &         ' ignored by MINIMIZE'
          call xerror(errmsg,80,1,21)
          nwarning = nwarning + 1
        END IF
        IF(md_respa) THEN
          errmsg=
     &         'Warning TIMESTEP(&INTEGRATOR) option is'/ /
     &         ' ignored by MINIMIZE.Timestep set to 1 fs'
          call xerror(errmsg,80,1,21)
          nwarning = nwarning + 1
        END IF
        md_respa=.FALSE.
        IF(thermos) THEN
          errmsg=
     &         'Warning THERMOS(&SIMULATION) option is'/ /
     &         ' ignored by MINIMIZE'
          call xerror(errmsg,80,1,21)
          nwarning = nwarning + 1
        END IF
        IF(cpress) THEN
          errmsg=
     &         'Warning (ISO)STRESS(&SIMULATION) option is'/ /
     &         ' ignored by MINIMIZE'
          call xerror(errmsg,80,1,21)
          nwarning = nwarning + 1
        END IF
        IF(maxstp .EQ. 0) THEN
          errmsg=
     &         'The number of minimization steps is zero. '/ /
     &         ' Use the TIME(&RUN) option to change steps.'
          call xerror(errmsg,80,1,20)
          nsevere = nsevere + 1 
        END IF
        IF(.NOT. stretch .OR. stretch_heavy) THEN
          errmsg=
     &         'Minimization cannot be run with constraints.'/ /
     &         ' Change STRETCH options.'
          call xerror(errmsg,80,1,20)
          nsevere = nsevere + 1 
        END IF
        IF(agbnp.and.solven) THEN
          errmsg= 
     &         'AGBNP mode: No solvent allowed  ' 
          call xerror(errmsg,80,1,21)
          nsevere = nsevere+1
        ENDIF  
        IF(agbnp) THEN 
          if(.not.(near0(alf-90.0).and.near0(bet-90.0).and.near0(gam
     &         -90.0)))THEN 
c     &        .not.near0(bet-90.0).and 
c     &        .not.near0(gam-90.0)) THEN  
            alf=90.d0 
            bet=90.d0 
            gam=90.d0 
            errmsg= 
     &           'AGBNP mode: cell angles reset to 90 degrees ' 
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
          END IF
        END IF
      END IF

      IF(frequencies .AND. mdsim) THEN
         errmsg=
     &'Calculation of harmonic frequencies is '/ /
     &'incompatible with MDSIM(&SIMULATION).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(frequencies .AND. md_respa) THEN
         errmsg=
     &'Calculation of harmonic frequencies is '/ /
     &'incompatible with &INTEGRATOR.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(anxrms .AND. nxrms .EQ. 0) THEN
         errmsg='inst_xrms(COMPARE(&PROPERTIES)): Cannot compute'
     & / /       ' without an output file.'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF

      IF(anxrms .AND. (.NOT. template)) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot compute X_rms without'
     &   / /     ' template file.'/ /' Provide TEMPLATE(&SETUP) cmd.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_str .AND. (.NOT. template)) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Running without' / /
     &        ' template file.'/ /' Provide TEMPLATE(&SETUP) cmd.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_ca .AND. avg_he) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot fit' / /
     &        ' structures to heavy and ca atoms at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_str .OR. avg_rms) THEN
         IF(.NOT. avg_ca .AND. .NOT. avg_he) THEN
            errmsg=
     &           'STRUCTURES(&PROPERTIES): Unknown atoms to which' / /
     &           ' fit structure. Use *averaged*.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
         END IF
      END IF
      IF(linser .AND. nflag(1) .NE. 0) THEN
         errmsg=
     &'CONTROL(&RUN) must be 0 when INSERT(&SET_UP) is'/ /
     &' specified' 
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

*=======================================================================
*---  Abort if control parameters are wrong (abort also atomic system)--
*=======================================================================
      
      IF(md_respa) THEN

         abort2=sction .OR. noneq
         IF(abort2) THEN
            errmsg=
     &'NON_EQUILIBRIUM(&SOLUTE) and MTS_RESPA(&INTEGRATOR)'/ /
     &'are mutually exclusive'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=rcutl-rcuth
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) radius_l > radius_h.'/ /
     &' Choose radius_l less or equal to radius_h'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=rcutm-rcutl
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) radius_m > radius_l.'/ /
     &'Choose radius_m less or equal to radius_l'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=km-kl
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) wrong k-ewald directive.'/ /
     &' Choose km less or equal to kl' 
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=kl-rkcut
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            write(kprint,10230)
10230       FORMAT(/ /  
     &'**ERROR: in MTS_RESPA(&INTEGRATOR) wrong k-ewald directive.' /
     &'- choose kl less or equal to rkcut in EWALD(&POTENTIAL)')
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(.NOT. md_respa) THEN
         IF(landersen) THEN
            errmsg=
     &'Directive ANDERSEN(&SIMULATION) can be used only'/ /
     &' when MTS_RESPA(&INTEGRATOR) is specified'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(nflag(1) .NE. 0 .AND. slv_create) THEN
         errmsg=
     &   'Cannot create solvent and read a restart files'
     & / /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(sgroup .AND. nflag(1) .NE. 0) THEN
         errmsg=
     &'Cannot apply symmetry and read a restart file'
     &/ /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      IF(sgroup .AND. tpgfil) THEN
         errmsg=
     &'Cannot read a binary TPGPRM file and apply symmetry'
     &/ /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(nplot_fragm.gt.0.and.nfragm.eq.0) THEN 
         errmsg=
     &   'No fragment defined in &READ_SOLUTE:' / /
     &   ' no fragment will be printed'
         call xerror(errmsg,80,1,21)
         nwarning= nwarning+1
      END IF   

      IF(slv_create .AND. (.NOT. tpgprm_read) ) THEN
         errmsg=
     &   'Cannot create solvent without '
     & / /'READ_PRM_ASCII and READ_TPG_ASCII in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_create .AND. nbun_slv .EQ. 0) THEN
         errmsg=
     &   'Cannot create solvent without '
     & / /'defining residues with JOIN SOLVENT in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. nbun_slv .EQ. 0) THEN
         errmsg=
     &   'Cannot add solvent without '
     & / /'defining residues with JOIN SOLVENT in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(pdb_read .AND. slv_create) THEN
         errmsg=
     &   'Cannot create solvent and read '
     & / /'from the PDB file defined in READ_PDB.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(pdb_read .AND. slt_create) THEN
         errmsg=
     &   'Cannot generate solvent and read '
     & / /'solute from the PDB file defined in READ_PDB.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(nflag(1).NE.0) THEN 
        if(.not.tpgfil) THEN 
          errmsg='Provide READ_TPGPRM directive when restarting'
     &         / / ' a simulation'
          call xerror(errmsg,80,1,20)
          nsevere=nsevere +1
        END IF
      END IF
      IF(.NOT. analys) THEN
         IF((.NOT. pdb_read .AND. .NOT. slt_create .AND. .NOT.
     &        slv_create).AND. nflag(1) .EQ. 0 .AND. .NOT. stoprun) THEN
            errmsg=
     &   'Cannot run without initial coordinates '
     & / /'from the PDB file defined in READ_PDB.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

      IF(kcoord_slt .EQ. 0 .AND. slt_create) THEN
         errmsg=
     &   'When creating solvent, solute coords. are read'
     & / /' by command COORDINATES in &SOLUTE.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(stoprun .AND. slv_create) THEN
         errmsg=
     &   'Cannot use STOP_RUN (&PARAMETERS) and'
     & / /' create solvent at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. sgroup) THEN
         errmsg=
     &   ' Command SGROUP (&SOLUTE) and ADD_UNITS (&SOLVENT) are'
     & / /' incompatible.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. slv_create) THEN
         errmsg=
     &   ' Cannot create solvent and use ADD_UNITS (&SOLVENT) '
     & / /'at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF


      IF(abmd .AND. DABS(spring) .LT. 1.0D-3) THEN
         errmsg=
     &   'In &POTENTIAL: The ABMD force constant is zero.'
     & / /' Hope this is ok..'
         call xerror(errmsg,80,1,21)
         nwarning= nwarning+1
      END IF

      IF(abmd .AND. diss_atoms) THEN
         errmsg=
     &   'In &POTENTIAL: In this implementation Cannot use'
     & / /'ABMD with subcommand (atoms).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(abmd .AND. (dissociate .AND. associate)) THEN
         errmsg=
     &   'In ABMD(&POTENTIAL): Cannot use'
     & / /' dissociate and associate at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(abmd) THEN
         ok=.TRUE.
         IF(fold .AND. dissociate) ok=.FALSE.
         IF(fold .AND. associate) ok=.FALSE.
         IF(fold .AND. abmd_tors) ok=.FALSE.
         IF(dissociate .AND. abmd_tors) ok=.FALSE.
         IF(associate .AND. abmd_tors) ok=.FALSE.
         IF(.NOT. ok) THEN
            errmsg=
     &           'In ABMD(&POTENTIAL): Cannot use'
     & / /' two ABMD options at the same time.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
         END IF
      END IF
      IF(coupl_mol .AND. coupl_grp) THEN
         errmsg=
     &        'In COUPLING(&SIMULATION): Can use only one type of '
     &     / /   'coupling.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF((coupl_atm.OR.coupl_grp) .AND. (.NOT.stretch)) THEN 
         errmsg=
     &        'Constraints on non X-H bonds are incompatible with '
     &     / /   'COUPLING group|atomic(&SIMULATION) option:'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
       END IF

*=======================================================================
*---  Error on properties ----------------------------------------------
*=======================================================================

      IF((.NOT. analys) .AND. near0(fconf) .AND. dmprnd) THEN
         errmsg=
     &   'Use TRAJECTORY(&INOUT) only with '
     &   / /' the environment &ANALYSIS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      IF(analys .AND. (.NOT. dmprnd)) THEN
         errmsg=
     &   'Cannot analyse a trajectory without'
     &   / /' reading from a file. Use TRAJECTORY(&INOUT).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
         
      IF((.NOT. analys) .AND. maxrun .EQ. 0 .AND. dmprnd ) THEN
         errmsg=
     &   'Cannot dump trajectory without the '
     &   / /'length of the run. Use MAXRUN(&RUN).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(analys .AND. gofr) THEN
         IF(gofr_neighbor .AND. update_anl .EQ. 0) THEN
            errmsg=
     &   'Cannot compute G of R with neighbor(GOFR(&ANALYSIS))'
     &   / /' without update frequency UPDATE(&ANALYSIS).'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      
      IF((hbonds_tot .OR. hbonds_res) .AND. nhbonds .EQ.
     &     0) THEN
         errmsg=
     &' Output file for hydrogen bonds computation is undefined. '
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(time_corr) THEN
         IF(not_time_corr) THEN
            errmsg=
     &'Can run TIME_CORRELATION(&PROPERTY) only without other '
     &           / /'properties.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

*===== Some option in properties are incompatible with md_respa =======

      IF(md_respa) THEN
         IF(vacf) THEN
            errmsg=
     &           ' VACF(&PROPERTIES) has no action while running'
     &           / /' a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(diffusion) THEN
            errmsg=
     &           ' DIFFUSION(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(fragm_dist) THEN
            errmsg=
     &           ' FRAGM_DIST(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(hbonds_tot .OR. hbonds_res) THEN
            errmsg=
     &           ' HBONDS(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(avg_rms) THEN
            errmsg=
     &           ' print rms in STRUCTURES(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
            
      END IF

      time=tsign*time

*=======================================================================
*---  Error counts is over; build abort message if errors are found 
*=======================================================================

c--   if errors are severe then 

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         idum1 = nsevere/10
         idum2 = nsevere - 10*idum1 
         idum1=idum1+1
         idum2=idum2+1
         if(idum2.eq.1) idum2=11
         nerr = integers(idum1:idum1)/ /integers(idum2:idum2)
         errmsg= nerr/ /' ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,2)
      ELSE IF(nsevere.gt.99) THEN 
         iret=1
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,2)
      END IF

c--   if warnings are given then 

      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0 
         idum1 = nwarning/10
         idum2 = nwarning - 10*idum1 
         idum1=idum1+1
         idum2=idum2+1
         if(idum2.eq.1) idum2=11
         nerr = integers(idum1:idum1)/ /integers(idum2:idum2)
         errmsg= nerr/ /' ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         iret=1
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,21)
         STOP
      END IF
      if(nsevere.eq.0.and.nwarning.eq.0) THEN
         write(kprint,'(5x,a)') 'Input OK!!'
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
