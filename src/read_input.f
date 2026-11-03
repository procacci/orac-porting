      SUBROUTINE read_input(sitew,iret,errmsg)

*======================================================================*
*                                                                      *
*                             Last update                              *
*                                                                      *
*     $Id: read_input.f,v 1.1.1.1 2009/10/19 07:48:42 signo Exp $
*                                                                      *
*     Written by Massimo Marchi CE Saclay                              *
*                                                                      *
*                                                                      *
*     EXTERNALS : Free format input package, MAXA0.                    *
*                                                                      *
*                                                                      *
************************************************************************


!======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar
      use rem
      use sge, only: read_sge
#ifdef _MPI_
      use orac_mpi !if used read_input knows everything 
#endif
#ifdef _OMP_
      use omp_integr
#else
      use omp_integr, only:omp_timing
#endif
#ifdef _OMP_
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#endif

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      INTEGER sitew,iret
      CHARACTER*80 errmsg

!-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER nword,ierr,naux
      LOGICAL lprint,lthr,lthr1,lthr2
      CHARACTER*8 fmt
      CHARACTER*22 err_open
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*15 err_fnf
      CHARACTER*80 line,strngs(40)
      CHARACTER*1 sep(2),comm(2)
      REAL*8  fmaxstp,fscale,fprint,fsave,fplot,fascii,fprop,fplot_fragm
     &     ,fplot_steer,fplottemp_stn0,fplot_stn0,fplot_alch0
     &     ,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms,fvi
     &     ,ftop_print,gofr_fprint,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf
     &     ,fnovaf,fdipole,favg,favg_xrms,ffragm_dist,fhbonds,frms,fgyr
     &     ,fabmd
      DATA sep/' ',','/comm/'#','!'/

!-------------------- COMMON VARIABLES ---------------------------------
#ifdef _OMP_
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#endif

      COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
      INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)

!==================== EXECUTABLE STATEMENTS ============================


!--------------  synatx errors messages---------------------------------
!   control error messages in all envs except PARAMETERS and ANALYSIS
  
      err_open = 'OPEN keyword not found'
      err_args(1) = 'Number of arguments must be at least '
      err_args(2) = 'Number of arguments must not exceed  '
      err_args(3) = 'Number of arguments must be '
      err_unr(1)='Unrecognized command  ---> '
      err_unr(2)='Unrecognized subcommand -> '
      err_unr(3)='Unrecognized keyword ----> '
      err_unr(4)='UNSUPPORTED  COMMAND ----> '
      err_end = '...or missing &END'  
      err_fnf = ' file not found' 

!------------- Set default values of some input parameters -------------
      CALL blkdta
      
      fmaxstp=0.0D0
      fscale=0.0D0
      fprint=0.0D0
      fsave=0.0D0
      fplot=0.0D0
      fascii=0.0D0
      fprop=0.0D0
      fplot_fragm=0.0D0
      fplot_steer=0.0D0
      fplottemp_stn0=0.0D0
      fplot_center=0.0D0
      frject=0.0D0
      fconf=0.0D0
      fmaxrun=0.0D0
      fupdte=0.0D0
      fxrms=0.0D0
      fvi=0.0D0
      ftop_print=0.0D0
      gofr_fprint=0.0D0
      gofr_favg=0.0D0
      gofr_fcomp=0.0D0
      fprtvaf=0.0D0
      ftotvaf=0.0D0
      fnovaf=0.0D0
      fdipole=0.0D0
      favg=0.0D0
      favg_xrms=0.0D0
      ffragm_dist=0.0D0
      fhbonds=0.0D0
      frms=0.0D0
      fgyr=0.0D0
      lthr=.false.
      lthr1=.false.
      lthr2=.false.
      lprint=.FALSE.
!     This is the copyright heading on s/o. Written by the root process 
#ifdef _MPI_
      IF (iproc.eq.0 ) THEN
        WRITE(*,*)
        WRITE(*,3)
        WRITE(*,*)
        WRITE(*,4)
        WRITE(*,*)
        WRITE(*,5)
        WRITE(*,*)
        WRITE(*,6)
        WRITE(*,*)
        WRITE(*,7)
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,10)
        WRITE(*,*)
        WRITE(*,11)
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,12)
        WRITE(*,*)
        WRITE(*,13)
        WRITE(*,*)
        WRITE(*,14)
        WRITE(*,*)
        WRITE(*,*)
      ENDIF
! print headings for the parallel run
#endif

      write(kprint,*)
      write(kprint,3)
      write(kprint,*)
      write(kprint,4)
      write(kprint,*)
      write(kprint,5)
      write(kprint,*)
      write(kprint,6)
      write(kprint,*)
      write(kprint,7)
      write(kprint,*)
      write(kprint,*)
      write(kprint,10)
      write(kprint,*)
      write(kprint,11)
      write(kprint,*)
      write(kprint,*)
      write(kprint,12)
      write(kprint,*)
      write(kprint,13)
      write(kprint,*)
      write(kprint,14)
      write(kprint,*)
      write(kprint,*)

      write(kprint,16)
101   continue
      line(79:80)='  '
      titld0(79:80)='  '

!----- Title card,$ 'blank'


c=======================================================================
c     Input parser starts here for file sys.mddata 
c=======================================================================


100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
#ifdef _OMP_ 
      IF(line(1:1) .EQ. '#')  THEN 
        if(line(1:2).EQ.'#&'.OR.line(1:3).EQ.'#&T'.AND.(.NOT.lthr)) THEN
          if(line(1:3).EQ.'#&T') omp_timing=.true.
          CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
          IF(iret.EQ.1) THEN 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            STOP
          END IF
          IF(strngs(2).EQ. 'NTHREADS') THEN
            IF(nword.eq.3) THEN 
              lthr=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr
!             print OMP headings      
            ELSE IF(nword.eq.5) THEN
              lthr=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr
              CALL fndfmt(1,strngs(5),fmt)
              READ(strngs(5),fmt,err=20) cache_line_size
            ELSE
              errmsg="Wrong threads num or cache_line_size"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'NT-LEVEL1') THEN 
            IF(nword.eq.3) THEN 
              lthr1=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr1
            ELSE IF(nword.eq.5) THEN
              lthr1=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr1
              CALL fndfmt(1,strngs(5),fmt)
              READ(strngs(5),fmt,err=20) cache_line_size1
            ELSE
              errmsg="Wrong threads num or cache_line_size"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'NT-LEVEL2') THEN 
            IF(nword.eq.3) THEN 
              lthr2=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr2
            ELSE IF(nword.eq.5) THEN
              lthr2=.true.
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) nthr2
              CALL fndfmt(1,strngs(5),fmt)
              READ(strngs(5),fmt,err=20) cache_line_size2
            ELSE
              errmsg="Wrong threads num or cache_line_size"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_INTEGR') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_integr=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_CNST') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_cnst=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_N0') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_n0=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_N1') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_n1=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_SHELL_H') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_shell_h=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_SHELL_L') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_shell_l=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_SHELL_M') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_shell_m=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_UPDATE') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_update=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE IF(strngs(2).EQ. 'LEVEL_FFTW') THEN 
            IF(nword.eq.3) THEN 
              CALL fndfmt(1,strngs(3),fmt)
              READ(strngs(3),fmt,err=20) naux
              level_fftw=naux
            ELSE
              errmsg="Expecting an integer in LEVEL args"
              CALL xerror(errmsg,80,1,30)
              STOP
            END IF
          ELSE
            errmsg='Expecting NTHREADS string  '
            CALL xerror(errmsg,80,1,30)
            STOP
          END IF
          go to 100
        ELSE
          go to 100
        END IF
      END IF
#else
      IF(line(1:1) .EQ. '#')  THEN 
        if(line(1:3).EQ.'#&T') omp_timing=.true.
        GOTO 100 
      END IF
#endif
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         go to 100
      END IF

c=====Environment RUN ==================================================
      IF(strngs(1).EQ. '&RUN' ) THEN
         CALL read_run(fprop,frject,fmaxrun,fmaxstp,fprint,err_args
     &        ,err_unr,err_end)
c=======================================================================

c=====Environment INOUT ================================================
      ELSE IF(strngs(1).EQ. '&INOUT' ) THEN
        CALL read_inout(fconf,fsave,fplot,fplot_fragm,fplot_steer
     &       ,fplottemp_stn0,fplot_stn0,fplot_alch0,fplot_center,fascii
     &       ,err_open,err_args,err_end,err_unr)
c=======================================================================

c=====Environment SIMULATION============================================
      ELSE IF(strngs(1).EQ. '&SIMULATION' ) THEN
         CALL read_simulation(fscale,err_args,err_unr,err_open,err_end)
c=======================================================================

c=====Environment POTENTIAL=============================================
      ELSE IF(strngs(1).EQ. '&POTENTIAL' ) THEN
         CALL read_potential(fupdte,fabmd,err_args,err_unr,err_end)
c=======================================================================

#ifdef _MPI_
c=====Environment REM===================================================
      ELSE IF(strngs(1).EQ. '&REM' ) THEN
         CALL read_rem(err_args,err_unr,err_end,err_fnf)
c=======================================================================
#endif

c=====Environment SGE===================================================
      ELSE IF(strngs(1).EQ. '&SGE' ) THEN
         CALL read_sge(err_args,err_unr,err_end,err_fnf,err_open)
c=======================================================================

c=====Environment META============================================
      ELSE IF(strngs(1).EQ. '&META' ) THEN
         CALL read_meta(err_args,err_unr,err_end)
c=======================================================================

c=====Environment SETUP=================================================
      ELSE IF(strngs(1).EQ. '&SETUP' ) THEN
         CALL read_setup(err_open,err_args,err_end,err_unr
     &        ,err_fnf)
c=======================================================================

c=====Environment INTEGRATOR============================================
      ELSE IF(strngs(1).EQ. '&INTEGRATOR' ) THEN
         CALL read_integrator(err_args,err_unr,err_end)
c=======================================================================

c=====Environment PARAMETERS============================================
      ELSE IF(strngs(1).EQ. '&PARAMETERS' ) THEN
         CALL read_parameters(iret,errmsg)
         IF(iret .NE. 0) THEN 
            CALL xerror(errmsg,80,1,2)
            STOP
         END IF
c=======================================================================

c=====Environment SOLUTE================================================
      ELSE IF(strngs(1).EQ. '&SOLUTE' ) THEN
         CALL read_solute(err_args,err_unr,err_end,err_fnf,err_open)
c=======================================================================

c=====Environment SOLVENT===============================================
      ELSE IF(strngs(1).EQ. '&SOLVENT' ) THEN
         CALL read_solvent(err_open,err_fnf,err_args,err_unr,err_end)
c=======================================================================

c=====Environment ANALYSIS==============================================
      ELSE IF(strngs(1).EQ. '&ANALYSIS' ) THEN
         CALL read_analysis(err_args,err_unr,err_end)

c=======================================================================

c=====Environment PROPERTIES ===========================================
      ELSE IF(strngs(1).EQ. '&PROPERTIES') THEN
         CALL read_properties(fxrms,favg,favg_xrms,fvi,gofr_fprint
     &        ,ftop_print,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf
     &        ,fnovaf,fdipole,ffragm_dist,fhbonds,frms,fgyr,err_open
     &        ,err_args,err_unr,err_end)
c=======================================================================

c=====Blank Line =======================================================
      ELSE IF(strngs(1).EQ. ' ' ) THEN
         CONTINUE

c=====Unrecognized Environment==========================================
      ELSE
         errmsg='Unrecognized ENVIRONMENT  ---> '// strngs(1)
         call xerror(errmsg,80,1,2) 
         STOP
      END IF
      
      GOTO 100

*-----------------------------------------------------------------------

600   CONTINUE

      write(kprint,17)
      CALL verify_input(fmaxstp,fscale,fprint,fsave,fplot,fascii
     &     ,fprop,fplot_fragm,fplot_steer,fplottemp_stn0,fplot_stn0
     &     ,fplot_alch0,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms
     &     ,fvi,favg,favg_xrms,ftop_print,gofr_fprint,gofr_favg
     &     ,gofr_fcomp,fprtvaf,ftotvaf,fnovaf,fdipole,ffragm_dist
     &     ,fhbonds,frms,fgyr,fabmd,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
102   CONTINUE 

#ifdef _OMP_
!      check that OpenMP execution choices are coherent 
      if((lthr1.OR.lthr2).and.(.not.lthr)) THEN 
        errmsg='NTHREADS (MAX) must be specified for OpenMP execution'
        call xerror(errmsg,80,1,2) 
        STOP
      end if
      if(nthr.gt.mcor3) THEN 
        errmsg='MAX OpenMP THREADS exceeded: change mcor3 in omp_integr'
        call xerror(errmsg,80,1,2) 
        STOP
      ENDIF
      if(nthr1.gt.mcor1) THEN 
        errmsg='MAX OpenMP THREADS exceeded: change mcor1 in omp_integr'
        call xerror(errmsg,80,1,2) 
        STOP
      ENDIF
      if(nthr2.gt.mcor2) THEN 
        errmsg='MAX OpenMP THREADS exceeded: change mcor2 in omp_integr'
        call xerror(errmsg,80,1,2) 
        STOP
      ENDIF
      if(lthr1.OR.lthr2) THEN 
        omp_dynamic=.true.
      ELSE
        omp_dynamic=.false.
        nthr1=nthr
        nthr2=nthr
      end if
      if(nthr.le.2) THEN 
        nthr1=nthr
        nthr2=nthr
        omp_dynamic=.false.
      END IF
      if(nthr1.ge.nthr.and.nthr2.ge.nthr) THEN 
        nthr1=nthr
        nthr2=nthr
        omp_dynamic=.false.
      END IF
!     set default behaviour for OMP executable
      if((.not.lthr) .AND. (.not.lthr1). AND. (.not.lthr2)) THEN 
!$OMP PARALLEL SHARED(nthr) 
        nthr=OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
        if(nthr.gt.mcor3) THEN 
          errmsg=
     &         'MAX OpenMP THREADS exceeded: change mcor3 in omp_integr'
          call xerror(errmsg,80,1,2) 
          STOP
        ENDIF
        if (nthr.gt.mcor1) THEN 
          nthr1=mcor1
        ELSE
          nthr1=nthr
        ENDIF
        if(nthr.gt.mcor2) THEN 
          nthr2=mcor2
        ELSE
          nthr2=nthr
        END IF
        omp_dynamic=.false.
      end if
      
#ifdef _PGI_
      omp_dynamic=.false.       ! no support for dyn thread in PGI
      nthr1=nthr
      nthr2=nthr
#endif
#ifndef _MPI_ 
      write(kprint,1010) nthr
     &     ,cache_line_size,nthr1,cache_line_size1,nthr2
     &     ,cache_line_size2,omp_dynamic,omp_timing
#else
!             print OMP-MPI headings      
      write(kprint,2010) ntrajectories
     &     ,nthr,cache_line_size,nthr1,cache_line_size1,nthr2
     &     ,cache_line_size2,omp_dynamic,omp_timing
#endif
      omp_dynamic=.false.       ! set false omp_dynamic in ANY case; no dynamic
                                ! threading *inside* parallel region
                                ! thread_num is fixed *before* entering in
                                ! the parallel region and it may adjusted
                                ! for each parallel region (e.g. force or mts_intra_n0/1)
#endif
      
#ifdef  _MPI_
!     print out headings on standard output
      call mpi_headings
#endif
      
      RETURN
      


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
1010            format
     &  (///10x,'====================================================' 
     &  /10x,'=                                                  =' 
     &  /10x,'=    O M P      L A Y E R    S T A R T E D         =' 
     &  /10x,'=                                                  =' 
     &  /10x,'=  3rd level MAX threads  ',i4,' cache_line',I4,'      =' 
     &  /10x,'=            1st level    ',i4,' cache_line',I4,'      ='   
     &  /10x,'=            2nd level    ',i4,' cache_line',I4,'      ='   
     &  /10x,'=            DYNAMIC THREADING ',L1,'                   ='   
     &  /10x,'=            OMP TIMING    ',L1,'                       =' 
     &  /10x,'===================================================='///)


2010            format
     &  (///10x,'====================================================' 
     &  /10x,'=                                                  =' 
     &  /10x,'=    O M P - M P I   O R A C   S T A R T E D       =' 
     &  /10x,'=                                                  =' 
     &  /10x,'=            NTRAJECTORIES ARE',i4,'                 =' 
     &  /10x,'=  3rd level MAX threads  ',i4,' cache_line',I4,'      =' 
     &  /10x,'=            1st level    ',i4,' cache_line',I4,'      ='   
     &  /10x,'=            2nd level    ',i4,' cache_line',I4,'      ='   
     &  /10x,'=            DYNAMIC THREADING ',L1,'                   ='   
     &  /10x,'=            OMP TIMING    ',L1,'                       =' 
     &  /10x,'===================================================='///)

1     FORMAT(80('='))
2     FORMAT('=',78x,'=')
3     FORMAT
     x    ('               oooooo       rrrrrrr        aaaaa    ',
     x     '      ccccccc'/
     x     '             oooooooooo     rrrrrrrr     aaaaaaaaa  ',
     x     '     ccccccccc'/
     x     '            ooo      ooo    rr    rrr    aaa   aaa  ',
     x     '    ccc'/
     x     '            oo        oo    rr    rr    aaa     aaa ',
     x     '    cc'/
     x     '            oo        oo    rr rrr      aaaaaaaaaaa ',
     x     '    cc'/
     x     '            ooo      ooo    rr  rrr     aaaaaaaaaaa ',
     x     '    ccc'/
     x     '             oooooooooo     rr    rr    aa       aa ',
     x     '     ccccccccc'/
     x     '               oooooo       rr    rrr   aa       aa ',
     x     '      ccccccc')
4     FORMAT('                              (Release 5)')
5     FORMAT(
     x'     "A Molecular Dynamics Program to Simulate',
     x' Complex Molecular Systems"')
6     FORMAT(
     x'            Copyright(C) 1989 - 2011  ',
     x' M. Marchi and P. Procacci')
7     FORMAT('                           All Right Reserved')
10    FORMAT(
     x'   ORAC is provided "as is" and without any',
     x' warranty express or implied.',/
     x'   The user assumes all risks of using ORAC.')
11    FORMAT(
     x'   The user may make copies of ORAC for his',
     x'/her own use, and modify those',/
     x'   copies. The user MAY distribute any copy of the',
     x' original or modified',/
     x'   source code to any users at any sites',
     x' under the General Public License')
12    FORMAT('   Department of Chemistry, University of Firenze,',
     x' Italy')
13    FORMAT('   http://www.chim.unifi.it/orac')
14    FORMAT('   Cite this program as:',/
     x'   S. Marsili, G. F. Signorini, R. Chelli, M. Marchi,',
     x' P. Procacci',/
     x'   J. Comput. Chem. 2010, 31, 1106-1116')
16    FORMAT(
     x     //'=========================== INPUT TO ',
     x     'THE RUN ==================================='/
     x     '=',78(' '),'=')
17    FORMAT('=',78(' '),'='/80('=')//)

#ifdef _OMP_
 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? Tab Character?? ' 
      CALL xerror(errmsg,80,1,2)
      STOP
#endif
      END 
