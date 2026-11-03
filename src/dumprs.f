      SUBROUTINE dumprs(kdump,restart_file,nstep,temp,ntap,ngrp,nprot
     &     ,xp0,yp0,zp0,vpx,vpy,vpz,xpg,ypg,zpg,xcm,ycm,zcm,vpx1,vpy1
     &     ,vpz1,vcx,vcy,vcz,gh,vh,vh1,nh,hoover,avg1,avg2,na1,na2,co,oc
     &     ,cpress,zz,zz1,lzz,lzz1,igr,krdf,igrn,save_all_restarts)

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
      use sge, only: dfi, estimate_df, df_av, ndf_av, dfi_inst,
     &     df_var_inst, cnt_upd, nstates, sge_update, sim_tempering

      IMPLICIT none

#ifdef _MPI_
      include 'mpif.h'
      integer :: ierr,para0,ntraj_eff
      
#endif
      
*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstep,ntap,kdump,ngrp,nprot
      INTEGER nh,na1,na2,igrn,krdf(*)
      LOGICAL hoover,cpress,igr,save_all_restarts
      REAL*8  xp0(*),yp0(*),zp0(*),vpx(*),vpy(*),vpz(*),xpg(*),ypg(*)
     &     ,zpg(*),vpx1(*),vpy1(*),vpz1(*),gh(*),vh(*),vh1(*),avg1(*)
     &     ,avg2(*),co(3,3),oc(3,3),zz(3,3),zz1(3,3),lzz,lzz1,temp,xcm(
     &     *),ycm(*),zcm(*),vcx(*),vcy(*),vcz(*)
      CHARACTER*80 restart_file

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k
      integer :: iproc, nproc

*==================== EXECUTABLE STATEMENTS ============================

#ifdef _MPI_
      call mpi_comm_size(mpi_comm_world, nproc, ierr)
      call mpi_comm_rank(mpi_comm_world, iproc, ierr)
#endif
      IF(save_all_restarts) THEN
#ifdef _MPI_
        if(rem_run) then
          ntraj_eff=nproc/nbatteries
          para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
          if(para0.eq.1) THEN 
            CLOSE(kdump)
            OPEN(unit=kdump,file=restart_file,form='FORMATTED')
          END IF
        else
          CLOSE(kdump)
          OPEN(unit=kdump,file=restart_file,form='FORMATTED')
        end if
#else
        CLOSE(kdump)
        OPEN(unit=kdump,file=restart_file,form='FORMATTED')
#endif 
      ELSE
        OPEN(unit=kdump,file=restart_file,form ='FORMATTED',status
     &       ='UNKNOWN')
        REWIND kdump
      END IF

#ifdef _MPI_
      if(rem_run) then
        if(.not.save_all_restarts) THEN 
          write(kdump,'(2i6)') para_index, nproc
        ELSE
!         do not write para_index and save only replica 1 if save all
!         restarts is true
          ntraj_eff=nproc/nbatteries
          para0=para_index-int((para_index-1)/ntraj_eff)*ntraj_eff
          if(para0.ne.1) RETURN
        END IF
      elseif(sim_tempering) then
#else
      nproc = 1
      iproc = 0
      if(sim_tempering) then
#endif
        write(kdump,'(4i6)') para_index, nstates, nproc, sge_update
        write(kdump,35) ( dfi(i), i=1, nstates-1 )
        write(kdump,'(10000i3)') ( estimate_df(i), i=1, nstates-1 )
        write(kdump,35) ( df_av(i), i=1, nstates-1 )
        write(kdump,35) ( ndf_av(i), i=1, nstates-1 )
        if( sge_update > 0 ) then
          do k = 1, sge_update
            write(kdump,35) ( dfi_inst(i,k), i=1, nstates-1 )
            write(kdump,35) ( df_var_inst(i,k), i=1, nstates-1 )
          enddo
          write(kdump,'(10000i6)') ( cnt_upd(i), i=1, nstates-1 )
        endif
      endif

35    format(10000g25.15)

      WRITE(kdump,*) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i=1,3)
      WRITE(kdump,*) nstep,temp
      WRITE(kdump,*) ntap,ngrp,nprot
      WRITE(kdump,*) (xp0(i),i=1,ntap),(yp0(i),i=1,ntap),
     &     (zp0(i),i=1,ntap),(vpx(i),i=1,ntap),(vpy(i),i=1,ntap)
     &     ,(vpz(i),i=1,ntap),(xpg(i),i=1,ngrp),(ypg(i),i=1,ngrp)
     &     ,(zpg(i),i=1,ngrp),(xcm(i),i=1,nprot),(ycm(i),i=1,nprot)
     &     ,(zcm(i),i=1,nprot),(vcx(i),i=1,nprot),(vcy(i),i=1,nprot)
     &     ,(vcz(i),i=1,nprot)
      WRITE(kdump,*) (vpx1(i),i=1,ntap),(vpy1(i),i=1,ntap),
     &     (vpz1(i),i=1,ntap)
      WRITE(kdump,*) na1,na2
      WRITE(kdump,*) (avg1(i),i=1,na1),(avg2(i),i=1,na2)
      WRITE(kdump,*) hoover,cpress

      IF(hoover) THEN
         WRITE(kdump,*) nh
         WRITE(kdump,*) (gh(i),i=1,nh),(vh(i),i=1,nh),(vh1(i),i=1,nh)
      END IF
      IF(cpress) THEN
        WRITE(kdump,*) ((zz(i,j),j=1,3),i=1,3),((zz1(i,j),j=1,3),i=1,3)
     &       ,lzz,lzz1
      END IF

      WRITE(kdump,*) igr
      IF(igr) THEN
         WRITE(kdump,*) igrn
         WRITE(kdump,*) (krdf(i),i=1,igrn)
      END IF
      CLOSE(kdump)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
