      PROGRAM main

************************************************************************
*                                                                      *
*                       ORAC of BLAKES 7                               *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*                              O R A C                                 *
*                                                                      *
*                           (5.2 Beta Version)                         *
*                                                                      *
*          "A Molecular Dynamics Program to Simulate Complex Systems"  *
*                                                                      *
*              Copyright(C) 1989 - 2011 by M.Marchi and P.Procacci     *
*                                                                      *
*                        All Right Reserved                            *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
* ORAC is provided "as is" and without any warranty express or implied.*
* The user assumes all risks of using ORAC.                            *
*                                                                      *
* The user may make copies of ORAC for his/her own use, and modify     *
* those copies. The user MAY NOT distribute any copy of the origina    *
* or modified source code to any users at any sites other than his/her *
* own.                                                                 *
*                                                                      *
*                                                                      *
*                                                                      *
*       Written by  Massimo Marchi and Piero Procacci                  *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar
      use rem
#ifdef _MPI_
      use orac_mpi
#endif
#ifdef _OMP_
#ifdef _BGQ_
      use omp_lib !omp layer 
#endif
      use omp_integr
#endif
      IMPLICIT none

*======================= DECLARATIONS ==================================

      INTEGER, allocatable :: mapnl(:,:),mapnlt(:,:)
   
      integer  mapnl_slv(nmap_slv)
      INTEGER  na,nb,nc,nd,ne,nf,ng,la
      PARAMETER(na=tsites, nb=tsitep, nc=elsiz, nd=tgroup, ne=hoov,
     x          nf=elsiz, ng=npm, la=1)
      REAL*8  xp0(nb),yp0(nb),zp0(nb),xpcm(ng),ypcm(ng),zpcm(ng),xpg(nd)
     &     ,ypg(nd),zpg(nd),gh(hoov)
      LOGICAL mapp(m1),mapp_slv(slvatm)

*-------------- DEFINITION OF THE SCRATCH COMMON BLOCK -----------------


*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER iret

#ifdef _OMP_
      INTEGER itask,istride,i,nmax
#ifndef _BGQ_
      include 'omp_lib.h'       ! OMP layer
#endif
#endif

*==================== EXECUTABLE STATEMENTS ============================
#ifdef _MPI_
! initialize the mpi environment
      call mpi_initialize
      
! create a different dir for each process 
! and put the current wdir in it 
      call mpi_mkdirs
      l_mpi=.true.
#endif
*=======================================================================
*----- Read form input files and initialiase topology arrays -----------
*=======================================================================
      
      allocate(mapnl(m8,1),stat=iret) ! will use the same memory for OMP version 
 
      CALL starta(mapnl,mapnl_slv,mapp,mapp_slv)


*=======================================================================
*----- Build the simulation box and, if needed, insert a macromolecule -
*----- in an equilibrated configuration of solvent molecules -----------
*=======================================================================

      CALL bldbox(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,iret,errmsg
     &     )

*----------------- If iret .EQ. 1 STOP ** ------------------------------
*                                 ====

      IF(iret .EQ. 2) CALL xerror(errmsg,80,1,20)
      IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      IF(.NOT. stoprun) THEN

*=======================================================================
*----- Verify some variables -------------------------------------------
*=======================================================================
        
        CALL verify_variables(iret,errmsg)
        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
        
*=======================================================================
*----- Run an MD simulation or analysis --------------------------------
*=======================================================================

#ifdef _OMP_ 
        nmax=4*m1
        call cache_line(nmax,nthr,array_omp,cache_line_size) ! initialize omp_array
        call cache_line(nmax,nthr1,array1_omp,cache_line_size1) ! initialize omp_array1
        call cache_line(nmax,nthr2,array2_omp,cache_line_size2) ! initialize omp_array2
        call OMP_SET_NUM_THREADS(nthr)  !set thread nukm to level 3
        m8t=1+m8/nthr
        allocate(mapnlt(m8t,nthr),stat=iret) 
!$OMP PARALLEL SHARED(istride) PRIVATE(itask)
        itask=1+OMP_GET_THREAD_NUM()
        istride=OMP_GET_NUM_THREADS()
        do i=1,m8t
          mapnlt(i,itask)=0
        end do
        call split_map(ntap,mapnl,grppt,ngrp,kprint,mapnlt(1,itask)
     &       ,itask,istride)
!$OMP END PARALLEL
        deallocate(mapnl,stat=iret) ! release memory for the big map 
        if ( iret /= 0 ) stop 'Deallocation failed for mapnl'

#endif
        IF(analys) THEN
          CALL drive_analysis(xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm
     &         ,ypcm,zpcm)
        ELSE IF(mdsim) THEN
#ifdef _OMP_
          CALL mtsmd(mapnlt,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm,ypcm,zpcm
     &         ,ntap)
#else
          CALL mtsmd(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm,ypcm,zpcm
     &         ,ntap)
#endif
        ELSE IF(minimize .OR. frequencies) THEN
#ifdef _OMP_
          CALL run_minimize(mapnlt,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm
     &         ,ypcm,zpcm)
#else
          CALL run_minimize(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm
     &         ,ypcm,zpcm)
#endif
        END IF
      END IF
      
      CALL close_files
#ifdef _MPI_
! finalize mpi run and print end-of-run message
      call mpi_close(iproc)
      WRITE(kprint,100)
#else
      WRITE(kprint,100)
#endif
      
*================= END OF EXECUTABLE STATEMENTS ========================
      
100   FORMAT(//
     &     '     *****************************************************',
     &     '*******************'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *                        PROGRAM COMPLETED           ',
     &     '                  *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *****************************************************',
     &     '*******************'/)
      
      END
      
      
