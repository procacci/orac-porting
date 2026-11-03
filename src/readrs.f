      SUBROUTINE readrs(kdump,nstep,temp,ntap,ngrp,nprot,xp0,yp0,zp0,vpx
     &     ,vpy,vpz,xpg,ypg,zpg,xcm,ycm,zcm,vpx1,vpy1,vpz1,vcx,vcy,vcz
     &     ,gh,vh,vh1,nh,hoover,avg1,avg2,na1,na2,co,oc,cpress,zz,zz1
     &     ,lzz,lzz1,igr,krdf,igrn,igrnmax,continue)

************************************************************************
*                                                                      *
*                                                                      *
*     READRS will read all data necessary to reastrt an MD run.        *
*                                                                      *
*     N1      :  Number of elements of AVG1 read by READRS.       (O)  *
*     N2      :  Number of elements of AVG2 read by READRS.       (O)  *
*     M1      :  Physical row dimension og IGR1 and IGR2.         (I)  *
*     IG1     :  Number of columns  of IGR1 read by READRS.       (O)  *
*     IG2     :  Number of columns  of IGR2 read by READRS.       (O)  *
*     MAXINT  :  Number of rows of IGR2 read by READRS.           (O)  *
*                                                                      *
*----- Last update 06/5/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS none.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use rem
      use sge, only: alloc2, dfi, estimate_df, df_av, ndf_av, dfi_inst,
     &     df_var_inst, cnt_upd, df, nstates, nstates_rs, sge_update,
     &     sge_update_rst, sim_tempering, zero_free_energy

      IMPLICIT none

#ifdef _MPI_
      include 'mpif.h'
      integer :: ierr
      logical :: exist
#endif
      integer :: pnrep, iproc, nproc

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstep,ntap,kdump,ngrp,m1,nprot
      INTEGER nh,na1,na2,igrn,igrnmax
      LOGICAL hoover,cpress,igr,krdf(*)
      REAL*8  xp0(*),yp0(*),zp0(*),vpx(*),vpy(*),vpz(*),xpg(*),ypg(*)
     &     ,zpg(*),vpx1(*),vpy1(*),vpz1(*),gh(*),vh(*),vh1(*),avg1(*)
     &     ,avg2(*),co(3,3),oc(3,3),zz(3,3),zz1(3,3),lzz,lzz1,temp,xcm(
     &     *),ycm(*),zcm(*),vcx(*),vcy(*),vcz(*)
      LOGICAL continue

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i, j, k, n1, n2
      LOGICAL igr1

*==================== EXECUTABLE STATEMENTS ============================
      
      REWIND kdump

#ifdef _MPI_
      call mpi_comm_size(mpi_comm_world, nproc, ierr)
      call mpi_comm_rank(mpi_comm_world, iproc, ierr)
#else
      nproc = 1
      iproc = 0
#endif

#ifdef _MPI_
      if(rem_run) then
! REM: read the index of the process (para_index) 
! and the total number of the processes (pnrep)
! para_index can be different from the rank of the actual rank of the
! process (iproc)
        read(kdump,*,end=100,err=100) para_index, pnrep
! the number of processes must not change
        if(pnrep.ne.nproc) then
          if(iproc.eq.0) then
            write(*,*)
            write(*,*)
     &           "*************************************************"
            write(*,*)
     &           "*               FATAL ERROR!!!                  *"
            write(*,*)
     &           "*                                               *"
            write(*,*)
     &           "* in REM the number of processors cannot change *"
            write(*,*)
     &           "*                                               *"
            write(*,*)
     &           "*************************************************"
          endif
          call mpi_finalize(ierr)
          stop
        endif
      elseif(sim_tempering) then
#else
      if(sim_tempering) then
#endif
! ST: read the index of the ensemble (para_index), the total number of ensembles (nstates_rs) and 
! the total number of the processes (pnrep)
        read(kdump,*,end=100,err=100) para_index, nstates_rs, pnrep,
     &       sge_update_rst
#ifdef _MPI_
        call mpi_barrier(mpi_comm_world,ierr)
#endif
        if(iproc == 0 ) then
          if(nstates_rs /= nstates) then
            write(*,*)
            write(*,*)
     &          "******************************************************"
            write(*,*)
     &          "*                   FATAL ERROR!!!                   *"
            write(*,*)
     &          "*                                                    *"
            write(*,*)
     &          "* In SGE simulations, the number of states from      *"
            write(*,*)
     &          "* the restart file does not correspond to the input. *"
            write(*,*)
     &          "* The program will stop.                             *"
            write(*,*)
     &          "*                                                    *"
            write(*,*)
     &          "******************************************************"
            write(*,*)
            write(*,*) 'Change SETUP(&SGE) accordingly.'
            write(*,*) 'The number of states must be: ', nstates_rs
#ifdef _MPI_
            call mpi_abort ( mpi_comm_world, ierr )
#else
            stop
#endif
          endif

          if(pnrep /= nproc) then
            write(*,*)
            write(*,*)
     & "*************************************************************"
            write(*,*)
     & "*                         WARNING!!!                        *"
            write(*,*)
     & "*                                                           *"
            write(*,*)
     & "* The number of replicas (processors) does not correspond   *"
            write(*,*)
     & "* to that used in the previous SGE run.                     *"
            write(*,*)
     & "* The simulation continues with the new number of replicas. *"
            write(*,*)
     & "*                                                           *"
            write(*,*)
     & "*************************************************************"
          endif
        endif

        call alloc2 ( iproc )

        read(kdump,*,end=100,err=100) ( dfi(i), i=1, nstates-1 )
        read(kdump,*,end=100,err=100) ( estimate_df(i), i=1, nstates-1 )
        read(kdump,*,end=100,err=100) ( df_av(i), i=1, nstates-1 )
        read(kdump,*,end=100,err=100) ( ndf_av(i), i=1, nstates-1 )

        if( sge_update_rst > 0 ) then
          if( sge_update > 0 ) then
            if( sge_update_rst > sge_update ) then
              do k = 1, sge_update
                read(kdump,*) ( dfi_inst(i,k), i=1, nstates-1 )
                read(kdump,*) ( df_var_inst(i,k), i=1, nstates-1 )
              enddo
              do k = sge_update+1, sge_update_rst
                read(kdump,*)
                read(kdump,*)
              enddo
            else
              do k = 1, sge_update_rst
                read(kdump,*) ( dfi_inst(i,k), i=1, nstates-1 )
                read(kdump,*) ( df_var_inst(i,k), i=1, nstates-1 )
              enddo
            endif
            read(kdump,*) ( cnt_upd(i), i=1, nstates-1 )
            do i = 1, nstates - 1
              if( cnt_upd(i) > sge_update ) cnt_upd(i) = 0
            enddo
          else
            do k = 1, 2*sge_update_rst
              read(kdump,*)
            enddo
            read(kdump,*)
          endif
        endif

        if( zero_free_energy ) then
          estimate_df(:) = 0
          df_av(:) = 0.d0
          ndf_av(:) = 0.d0
          if( sge_update > 0 ) then
            dfi_inst(:,:) = 0.d0
            df_var_inst(:,:) = 0.d0
            cnt_upd(:) = 0.d0
          endif
        endif

        df(:) = dfi(:)
        do i = 1, nstates-1
          if( ndf_av(i) /= 0.d0 ) then
            dfi(i) = df_av(i) / ndf_av(i)
            df(i) = dfi(i)
          endif
        enddo

      endif


      READ(kdump,*,END=100) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i
     &     =1,3)
      READ(kdump,*,END=400,ERR=400) nstep,temp
      GOTO 500
400   CONTINUE
      WRITE(6,*)
      WRITE(6,'(a)') ' An old restart file was found. I hope it is '//
     &               ' Ok.'
      WRITE(6,*)
500   CONTINUE

      READ(kdump,*,END=100,ERR=100) ntap,ngrp,nprot
      READ(kdump,*,END=100,ERR=100) (xp0(i),i=1,ntap),(yp0(i),i=1,ntap)
     &     ,(zp0(i),i=1,ntap),(vpx(i),i=1,ntap),(vpy(i),i=1,ntap)
     &     ,(vpz(i),i=1,ntap),(xpg(i),i=1,ngrp),(ypg(i),i=1,ngrp)
     &     ,(zpg(i),i=1,ngrp),(xcm(i),i=1,nprot),(ycm(i),i=1,nprot)
     &     ,(zcm(i),i=1,nprot),(vcx(i),i=1,nprot),(vcy(i),i=1,nprot)
     &     ,(vcz(i),i=1,nprot)
      READ(kdump,*,END=100,ERR=100) (vpx1(i),i=1,ntap),(vpy1(i),i=1
     &        ,ntap),(vpz1(i),i=1,ntap)

      READ(kdump,*,END=100,ERR=100) n1,n2

      IF(na1 .LT. n1 .OR. na2 .LT. n2) THEN
         WRITE(6,'(a)') 'In READRS: Dimension of the averages array'
     &        //' insufficient. Abort.'
         STOP
      END IF
      READ(kdump,*,END=100,ERR=100) (avg1(i),i=1,n1),(avg2(i),i=1,n2)

      READ(kdump,*,ERR=200,END=200) hoover,cpress
      IF(hoover) THEN
         READ(kdump,*,END=100,ERR=100) nh
         READ(kdump,*,END=100,ERR=100) (gh(i),i=1,nh),(vh(i),i=1,nh)
     &        ,(vh1(i),i=1,nh)
      END IF
      IF(cpress) THEN
        READ(kdump,*,END=100,ERR=100) ((zz(i,j),j=1,3),i=1,3),((zz1(i
     &       ,j),j=1,3),i=1,3),lzz,lzz1
      END IF

      READ(kdump,*,END=100,ERR=100) igr1
      IF(igr .AND. igr1) THEN
         READ(kdump,*,END=100,ERR=100) igrn
         IF(igrn .LE. igrnmax) THEN
            READ(kdump,*,END=100,ERR=100) (krdf(i),i=1,igrn)
         ELSE
            GOTO 300
         END IF
      END IF

      GOTO 200
100   CONTINUE
      IF(continue) THEN
         WRITE(6,*)
         WRITE(6,'(a)') ' End of file encountered. Wrong restart'//
     & ' file or file has been corrupted, but program will not stop.'
         WRITE(6,'(a)') '     Hope it is ok. '
         WRITE(6,*)
      ELSE
         WRITE(6,*)
         WRITE(6,'(a)') ' End of file encountered. Wrong restart'//
     &   ' file or file has been corrupted. Abort.'
         WRITE(6,*)
         STOP
      END IF

200   CONTINUE

      RETURN

300   CONTINUE
      WRITE(6,*)
      WRITE(6,'(a)') ' The G(r) is read in incorrectly. I hope it is'//
     &               ' ok. Run continues'
      WRITE(6,*)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
