      SUBROUTINE drive_analysis(xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm
     &     ,ypcm,zpcm)

************************************************************************
*   Time-stamp: <98/03/20 16:21:38 marchi>                             *
*                                                                      *
*     drive_analysis analize a trajectory file written by mtsmd        *
*     In addition to that file also a binary topology file must        *
*     be provided                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  2 1997 -                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*  drive_analysis externals:   	       	       	       	       	       *
*       appbou asng_xrmsw calc_avg2_str calc_avg_str                   *
*       calc_avg_xrms calc_gofr calc_xrms change_frame change_step_sim *
*       check_length_sim check_read_columns check_topology 	       *
*       comp_cell_center comp_dip comp_molmass           	       *
*       comp_rmsq coordinate_spline coordinate_spline_init             *
*       daxpy dcffti dcopy dscal                find_length_run fndgrp *
*       get_spectra_vacf get_type_slv get_velocities inicmp matinv     *
*       mts_plot_fragm plotc plotd print_title_analysis		       *
*       prtat prtba prtbnd prtit prtpt prtsq			       *
*       read_confc_columns read_confc_rows set_hbonds_masks timer      *
*       time_correlation tr_inbox update write_bends write_bonds       *
*       write_diffusion write_fragm_dist write_gofrp write_gofrw       *
*       write_gyr write_hbonds write_rms write_tors write_vacf_phi     *
*       write_xrms write_xrms_atm xerror zero0 zeroa zero_gofr         *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar
      use giofar
      use fourier
      use spme

      IMPLICIT none
      
      INTEGER ma,mb,mh,ms,mf,memp,t1,mbs,numpr
      PARAMETER (ma=tsites,mb=tsitep,mh=elsiz,ms=secsiz,memp=1
     &     ,mf=mb*4,t1=types*(types+1)/2)

      PARAMETER (mbs=2*mb, numpr=npm)

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),eta(*),xpcm(*)
     &     ,ypcm(*),zpcm(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      include 'lc_list.h'
      INCLUDE 'iobuffer.h'

*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      EXTERNAL  near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg

      INTEGER    l2
      INTEGER    nato_slt,iter_avg,iter_avg2
      REAL*8  fact,dtvi

      INTEGER ntot_fragm,fragm_1,fragm_2
      INTEGER i,j,nstep,iret,istep
      LOGICAL end

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER nnlpp0(mpp),nnlppf(sitslu+1)

      INTEGER mapdn(2,mf),nmapdn(mf)
      REAL*8  vpx(mb),vpy(mb),vpz(mb),vpx1(mb),vpy1(mb),vpz1(mb)
     &     ,vcax(numpr),vcay(numpr),vcaz(numpr),vi(mb),aux 
      REAL*8  vh1(hoov),vco(3,3),qt0(4)
      LOGICAL a_mask(m1),d_mask(m1)

*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external

      COMMON /dynam/ w1,w2,vpx,vpy,vpz,vpx1,vpy1,vpz1,vh1
     &     ,vcax,vcay,vcaz,vco,mapdn,nmapdn,nnlpp0,a_mask,d_mask
     &     ,xpc,ypc,zpc,vxpc,vypc,vzpc,wsave1,spline_x,spline_y
     &     ,vacf_data,rms_disp,tot_rms_disp,xpb,ypb,zpb,xpcc,ypcc,zpcc

      
      INTEGER krdf(maxint*g1),type_slv(slvatm)
      INTEGER ntype_slv,offset_slv,nbetab_slv
      PARAMETER(nbetab_slv=90)
      CHARACTER*1 betab_slv(nbetab_slv)
      INTEGER itype_slv(nbetab_slv)
      REAL*8  xpo(mb),ypo(mb),zpo(mb)
      REAL*4  xau(mb),yau(mb),zau(mb)
      REAL*8  xpa(mb),ypa(mb),zpa(mb),xpga(m11),ypga(m11),zpga(m11)
     &     ,xpcma(npm),ypcma(npm),zpcma(npm),fstep,fstep1,fstep2
      REAL*8  wca(m1),whe(m1),wbc(m1),errca(npm),errhe(npm
     &     ),errbc(npm),erral(npm),drpca(m1),drpbc(m1),drphe(m1)
     &     ,drpal(m1),xp_avg(m1),yp_avg(m1),zp_avg(m1),xp_avg2(m1)
     &     ,yp_avg2(m1),zp_avg2(m1),tmass(numpr),vfcp,tfcp,gcpu
      REAL*8 dips(3),vol_gofr,sum_volume,elapse

      INTEGER offset,nnstep,start_time,end_time,length_run,length_tot
     &     ,length_fft
     &     ,iatom,iatom0,niatom,naux,i_old,i_start,i_end,typei
     &     ,nato1,nato2,idum

      COMMON /rag2/ xpo,ypo,zpo,xau,yau,zau,xpa,ypa,zpa,xpga,ypga,zpga
     &     ,xpcma,ypcma,zpcma,wca,whe,wbc,errca,errhe,errbc,erral,drpca
     &     ,drpbc,drphe,drpal,xp_avg,yp_avg,zp_avg,xp_avg2,yp_avg2
     &     ,zp_avg2,tmass

*==================== EXECUTABLE STATEMENTS ============================

*===  Check if the dimension of the work array are sufficient 

      IF(mb.LT.ntap) THEN
         errmsg=' While running Analysis: PARAMETER MB dimensions the'
     &//' work arrays in sufficiently. Abort. '
         CALL xerror(errmsg,80,1,2)
      END IF
    
*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

      rspoff=rspcut
      rspcut=0.0D0

*===  set few variable to zero

      iter_avg=0
      iter_avg2=0

      sum_volume=0.0D0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv


      CALL comp_molmass(nprot,protl,mass,tmass)

      IF(hbonds_tot .OR. hbonds_res) THEN
         acut_hb=DCOS(acut_hb*pi/180.0D0)
         a2cut_hb=DCOS(a2cut_hb*pi/180.0D0)
      END IF

      IF(anxrms .OR. gofr .OR. avg_str .OR. (diffusion .AND. slt_exist))
     &     CALL asng_xrmsw(ss_point,m1+1,wca,whe,wbc,beta,mback,nbone)
      
      IF(anxrms_cell) THEN
         anprot=.TRUE.
         annpro=annpro+1
         anpoint(1,annpro)=1
         anpoint(2,annpro)=ntap
      END IF

      IF(gofr) CALL get_type_slv(nato_slt,nato_slv,beta,betab_slv
     &     ,nbetab_slv,ntype_slv,type_slv,itype_slv,offset_slv
     &     ,types_gofr,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

      IF(avg_str .OR. avg_rms .OR. (diffusion .AND. slt_exist)) CALL
     &     zeroa(xp_avg,yp_avg,zp_avg,ntap,1)
      IF(avg_rms) CALL zeroa(xp_avg2,yp_avg2,zp_avg2,ntap,1)

*=======================================================================
*---- Initialize G of Rs -----------------------------------------------
*=======================================================================

      IF(gofr) THEN
         CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
      END IF

      IF(nplot_fragm .GT. 0) THEN
         do i=1,nfragm
            IF ( fragm(1,i).gt.ntap.or.fragm(2,i).gt.ntap) THEN 
               WRITE(kprint,10260) i
10260          FORMAT(//  
     &              '**ERROR: fragments limits for',i5,' DEF_FRAGM',
     &              '(&READ_SOLUTE) exceeds protein atoms' )
               STOP
            END IF 
         end do
         ntot_fragm=0
         do i=1,nfragm
            fragm_1=fragm(1,i)
            fragm_2=fragm(2,i)
            ntot_fragm = ntot_fragm + fragm_2-fragm_1+1 
         end do    
      END IF

      if(ncofactor.gt.0.and.nvi.gt.0) then
c--      set field cut-off to outer shell radius and alphaf 
c--      to current value if non set    
         if(alphaf.lt.0.d0.or.cut_field.lt.0.d0) then
            alphaf=alphal
            cut_field=rcuth+rtolh
         end if
         if(ncofactor.ne.12345678) THEN 
           CALL field_init(ncofactor,nsegs,kprint,naux) 
           if(naux.eq.0) THEN
             ncofactor = 0
             write(kprint,10977) 
           endif
         ELSE
           ncofactor=ntap
           do i=1,ntap 
             nsegs(1,i)=i
             nsegs(2,i)=i
           end do
         ENDIF
         dtvi = time/float(lrespa*mrespa)
         cut_field=cut_field**2
         if(lform_field) THEN 
           write(kvi,10103) dtvi,ncofactor
10103      FORMAT(f20.5,I10)
           if(ncofactor.ne.ntap) THEN 
             do i=1,ncofactor
               write(kvi,10104) nsegs(1,i),nsegs(2,i)
10104          FORMAT(2I10) 
             end do
           END IF
         END IF
      end if

*=======================================================================
*---  set fake neighbor lists
*=======================================================================

      nnlppf(1)=ngrp
      DO j=1,ngrp
         nnlppf(1+j)=j
      end do

*=======================================================================
*---  Check if there are no error in the list of bonds, bendings -------
*---  and torsions -----------------------------------------------------
*=======================================================================

      IF(prttopl) THEN
         CALL check_topology('B',top_bonds,lbond,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('D',top_bendings,lbend,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('P',top_ptors,ltors,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('I',top_itors,litor,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*----- Find out what was the timestep ----------------------------------
*=======================================================================

      IF(dmprnd) THEN
         CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap,fstep1,1
     &        ,end,divide_records,atom_record)
         IF(end) THEN
            errmsg=' Trajectory file seems to be empty.'
            CALL xerror(errmsg,80,1,2)
         END IF
         CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap,fstep2,2
     &        ,end,divide_records,atom_record)
         IF(end) THEN
            errmsg=' Trajectory file contains only one configuration.'
            CALL xerror(errmsg,80,1,1)
            fstep=fstep1
            stop_anl=1
         ELSE
            fstep=fstep2-fstep1

*=======================================================================
*----- Find out the length of the run ----------------------------------
*=======================================================================

            naux=stop_anl
            CALL find_length_run(stop_anl)
            WRITE(kprint,2000) stop_anl
            IF(naux .LT. stop_anl) stop_anl=naux
         END IF
      ELSE
         errmsg=
     &' Trajectory files must be direct access only '
         CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*---  Write a Banner Page                                            ---
*=======================================================================

      WRITE(kprint,'(//)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1100)
      WRITE(kprint,1100)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(//)')

      IF(prttpg) THEN
         WRITE(kprint,1000)
         IF(prtseq) CALL prtsq(nbun,mend,prsymb)
         IF(prtatl) CALL prtat(ss_index,ntap,beta,betb,nres,m1,prsymb
     &        ,chrge,mass)
         IF(prtbndl) CALL prtbnd(beta,nres(1,1),lbnd,lbond,potbo
     &        ,m9)
*         IF(prtcnl) CALL prtcn(dssp,lcnstr,lconstr,beta,nres(1,1))
         IF(prtbal) CALL prtba(beta,nres(1,1),lbndg,lbend,potbe,m2)
         IF(prtptl) CALL prtpt(beta,nres(1,1),ltor,ltors,potto,m3)
         IF(prtitl) CALL prtit(beta,nres(1,1),litr,litor,potit,m4)
         WRITE(kprint,'(////77(''=''))')
      END IF

*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================

      CALL print_title_analysis(fstep)

*=======================================================================
*----- Decide where to start to read -----------------------------------
*=======================================================================

      nstep=0
      IF(start_anl .NE. 0) nstep=start_anl-1

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

      IF(time_corr) THEN

         IF(vacf) THEN
            WRITE(kprint,2500)
         END IF
         IF(diffusion) THEN
            WRITE(kprint,2501)
         END IF

*=======================================================================
*--- Check the length of the simulation --------------------------------
*=======================================================================

         CALL check_length_sim(start_anl,stop_anl,buffer_time,buffer_fft
     &        ,start_time,end_time,length_run,length_tot,length_fft,iret
     &        ,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--- Check if buffer is sufficiently large for simulation parameters ---
*=======================================================================

         CALL check_read_columns(nbuffer,length_tot,atom_record,iret
     &        ,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--- Initialize cubic spline arrays to compute velocities --------------
*=======================================================================

         CALL coordinate_spline_init(spline_x,fstep,length_tot)

*=======================================================================
*--- Reinitialize length of simulation and fft's if necessary ----------
*=======================================================================

         IF(vacf .AND. divide_spline .GT. 1) THEN
            CALL change_step_sim(buffer_time,buffer_fft,fstep
     &           ,length_tot,length_fft,divide_spline,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

*=======================================================================
*--- Initialize FFT routines -------------------------------------------
*=======================================================================

         CALL dcffti(length_fft,wsave1)

*=======================================================================
*--- Zero correlation arrays -------------------------------------------
*=======================================================================

         DO j=1,2
            CALL zero0(vacf_data(0,j),length_tot)
            CALL zero0(tot_rms_disp(0,j),length_tot)
         END DO

         IF(diffusion) THEN
            CALL zero0(xpcc,length_tot)
            CALL zero0(ypcc,length_tot)
            CALL zero0(zpcc,length_tot)

            WRITE(kprint,5000)
            nnstep=0
            DO istep=start_time,end_time
               nnstep=nnstep+1

*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
            
               CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap
     &              ,fstep1,istep,end,divide_records,atom_record)

               CALL comp_cell_center(ntap,xpcc(nnstep),ypcc(nnstep)
     &              ,zpcc(nnstep),xp0,yp0,zp0)

            END DO
         END IF

*=======================================================================
*----  Set atom number of solute and solvent to zero and decide --------
*----  for which atom to compute correlation ---------------------------
*=======================================================================

         nato1=0
         nato2=0

         IF(corr_atoms(1) .EQ. 0) THEN
            niatom=ntap
         ELSE
            niatom=corr_atoms(1)
         END IF

*=======================================================================
*---- Beginning of time correlation loop -------------------------------
*=======================================================================

         i_old=0
         DO iatom0=1,niatom
            IF(corr_atoms(1) .EQ. 0) THEN
               iatom=iatom0
            ELSE
               iatom=corr_atoms(1+iatom0)
            END IF
            typei=ss_index(iatom)

*=======================================================================
*---- Read coordinates of each atom as a function of time --------------
*=======================================================================

            CALL read_confc_columns(xpc,ypc,zpc,xpb,ypb,zpb,ntap,iatom
     &           ,length_run,start_time,end_time,divide_records
     &           ,atom_record,iret,errmsg)
            IF(vacf) THEN

*=======================================================================
*---- Compute velocity autocorrelation ---------------------------------
*=======================================================================

*---- Do x component

               CALL coordinate_spline(spline_x,spline_y,xpc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vxpc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vxpc,vxpc,wsave1
     &              ,vacf_data(0,typei),fstep,w1,w2)
               
*---- Do y component
               
               CALL coordinate_spline(spline_x,spline_y,ypc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vypc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vypc,vypc,wsave1
     &              ,vacf_data(0,typei),fstep,w1,w2)
               
*---- Do z component
               
               CALL coordinate_spline(spline_x,spline_y,zpc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vzpc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vzpc,vzpc,wsave1
     &              ,vacf_data(0,typei),fstep,w1,w2)
               
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

            END IF
            IF(diffusion) THEN
               IF(typei .EQ. 1) THEN
                  nato1=nato1+1
               ELSE IF(typei .EQ. 2) THEN
                  nato2=nato2+1
               END IF

*=======================================================================
*---- Compute rms displacement -----------------------------------------
*=======================================================================

               CALL zero0(rms_disp(0,typei),length_tot)

*---- Do x component
               
               CALL daxpy(length_tot,1.0D0,xpcc,1,xpc,1)
               CALL dcopy(length_tot,xpc,1,vxpc,1)
               CALL time_correlation(length_fft,vxpc,vxpc,wsave1
     &              ,rms_disp(0,typei),fstep,w1,w2)
               
*---- Do y component
               
               CALL daxpy(length_tot,1.0D0,ypcc,1,ypc,1)
               CALL dcopy(length_tot,ypc,1,vypc,1)
               CALL time_correlation(length_fft,vypc,vypc,wsave1
     &              ,rms_disp(0,typei),fstep,w1,w2)
               
*---- Do z component
               
               CALL daxpy(length_tot,1.0D0,zpcc,1,zpc,1)
               CALL dcopy(length_tot,zpc,1,vzpc,1)
               CALL time_correlation(length_fft,vzpc,vzpc,wsave1
     &              ,rms_disp(0,typei),fstep,w1,w2)
               
*--- Compute the displacement for the current atom and accumulate
               
               CALL comp_rmsq(length_tot,rms_disp(0,typei),xpc,ypc
     &              ,zpc,tot_rms_disp(0,typei),w1,w2)
            END IF
            naux=MOD(iatom,atom_record)
            IF(naux .EQ. 0) THEN
               naux=iatom/atom_record
            ELSE
               naux=iatom/atom_record+1
            END IF
            i_start=(naux-1)*atom_record+1
            i_end=naux*atom_record
            IF(i_start .NE. i_old) THEN
               IF(i_end .GT. ntap) i_end=ntap
               WRITE(kprint,3000) i_start,i_end
            END IF
            i_old=i_start
         END DO

         IF(vacf) THEN
            CALL get_spectra_vacf(buffer_time,length_tot,length_tot/2
     &           ,vacf_data,vacf_spectra,wsave1)
            WRITE(kprint,3500)
            CALL write_vacf_phi(buffer_time,length_tot,fstep,vacf_data
     &           ,vacf_spectra)
         END IF
         IF(diffusion) THEN
            WRITE(kprint,4500)
            CALL write_diffusion(buffer_time,length_tot,fstep
     &           ,tot_rms_disp,nato1,nato2)
         END IF
      ELSE
         nnstep=0
         DO WHILE(.NOT. end .AND. nstep .LT. stop_anl)
            nstep=nstep+1
            nnstep=nnstep+1
            
*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
            
            CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap,fstep
     &           ,nstep,end,divide_records,atom_record)
            IF(.NOT. end) THEN

*=======================================================================
*----- Write information about the processed trajectory ----------------
*=======================================================================
               
               IF(gofr . AND. update_anl .NE. 0) THEN
                  IF(MOD(nnstep-1,update_anl) .NE. 0) THEN
                     WRITE(kprint,90000) fstep
                  END IF
               ELSE 
                  WRITE(kprint,90000) fstep
               END IF

*=======================================================================
*----- Check if all coordinates are non zeroes -------------------------
*=======================================================================
               
               CALL check_zero_coord(xp0,yp0,zp0,ntap,iret,errmsg)
               IF(iret .NE. 0) call xerror(errmsg,80,1,21)
               
*=======================================================================
*--- Compute the volume of the system ----------------------------------
*=======================================================================
               
               CALL matinv(3,3,co,oc,volume)
               volume=volume*boxl**3
               sum_volume=sum_volume+volume
               
*=======================================================================
*-------- Calculate group position  ------------------------------------
*=======================================================================
               
               CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,ngrp,grppt)
               
*=======================================================================
*-------- Calculate solute center of mass ------------------------------
*=======================================================================
               
               CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &              ,nprot,protl)
               
*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================
               
               CALL fndgrp(nprot,protl,atomp)
               
*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================
               
               CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
               CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga
     &              ,zpga)
               CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &              ,ypcma,zpcma)
               IF(gofr) THEN
                  IF(update_anl .NE. 0) THEN
                     IF(MOD(nnstep-1,update_anl) .EQ. 0) THEN
                        CALL update(co,xpga,ypga,zpga,xpa,ypa,zpa
     &                       ,rspcut,rspoff,ngrp,nnlpp0,mpp)
                        WRITE(kprint,90000) fstep
                     END IF
                  END IF
               END IF
               
*=======================================================================
*---- Compute instantaneous X-rms --------------------------------------
*=======================================================================
               
               IF(anxrms) THEN
                  CALL calc_xrms(anxca,anxbc,anxhe,anxal,anprot,annpro
     &                 ,anpoint,protl,wca,whe,wbc,xp0,yp0,zp0,xpt0,ypt0
     &                 ,zpt0,nato_slt,errca,errhe,errbc,erral,drpca
     &                 ,drpbc,drphe,drpal)

                  IF(MOD(nstep,nxrms) .EQ. 0) THEN
                     IF(anprot) THEN
                        WRITE(kxrms,50000) fstep
                        IF(anxca) CALL write_xrms(kxrms,annpro,'CA'
     &                       ,errca)
                        IF(anxbc) CALL write_xrms(kxrms,annpro,'BC'
     &                       ,errbc)
                        IF(anxhe) CALL write_xrms(kxrms,annpro,'HE'
     &                       ,errhe)
                        IF(anxal) CALL write_xrms(kxrms,annpro,'AL'
     &                       ,erral)
                     ELSE
                        WRITE(kxrms,50000) fstep
                        IF(anxca) CALL write_xrms(kxrms,nprot,'CA',errca
     &                       )
                        IF(anxbc) CALL write_xrms(kxrms,nprot,'BC',errbc
     &                       )
                        IF(anxhe) CALL write_xrms(kxrms,nprot,'HE',errhe
     &                       )
                        IF(anxal) CALL write_xrms(kxrms,nprot,'AL',erral
     &                       )
                     END IF
                     REWIND kxrms_atm
                     IF(anxca) CALL write_xrms_atm(kxrms_atm,ntap,'CA'
     &                    ,drpca,nnstep,fstep,ngrp,grppt,nres(1,1))
                     IF(anxbc) CALL write_xrms_atm(kxrms_atm,ntap,'BC'
     &                    ,drpbc,nnstep,fstep,ngrp,grppt,nres(1,1))
                     IF(anxhe) CALL write_xrms_atm(kxrms_atm,ntap,'HE'
     &                    ,drphe,nnstep,fstep,ngrp,grppt,nres(1,1))
                     IF(anxal) CALL write_xrms_atm(kxrms_atm,ntap,'AL'
     &                    ,drpal,nnstep,fstep,ngrp,grppt,nres(1,1))
                  END IF
               END IF
               
               
*=======================================================================
*---- Compute averaged structure ---------------------------------------
*=======================================================================
               
               IF(avg_str .OR. avg_rms) THEN
                  IF(avg_ca) CALL calc_avg_str(anxrms_cell,protl,wca
     &                 ,xpt0,ypt0,zpt0,xp_avg,yp_avg,zp_avg,qt0,xp0,yp0
     &                 ,zp0,nato_slt,iter_avg)
                  IF(avg_he) CALL calc_avg_str(anxrms_cell,protl,whe
     &                 ,xpt0,ypt0,zpt0,xp_avg,yp_avg,zp_avg,qt0,xp0,yp0
     &                 ,zp0,nato_slt,iter_avg)
               END IF
               
               IF(navg_str .NE. 0) THEN
                  IF(MOD(nstep,navg_str) .EQ.0) THEN
                     CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                     CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                     CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                     fact=1.0D0/DFLOAT(iter_avg)
                     CALL dscal(nato_slt,fact,xpo,1)
                     CALL dscal(nato_slt,fact,ypo,1)
                     CALL dscal(nato_slt,fact,zpo,1)
                     IF(avg_ca) WRITE(kavg,80000) 
                     IF(avg_he) WRITE(kavg,80100) 
                     CALL plotd(fstep,kavg,beta,xpt0,ypt0,zpt0,xpo,ypo
     &                    ,zpo,nato_slt,nres,m1,prsymb)
                  END IF
               END IF
               
               IF(navg_str_xrms .NE. 0) THEN
                  IF(MOD(nstep,navg_str_xrms) .EQ.0) THEN
                     CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                     CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                     CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                     fact=1.0D0/DFLOAT(iter_avg)
                     CALL dscal(nato_slt,fact,xpo,1)
                     CALL dscal(nato_slt,fact,ypo,1)
                     CALL dscal(nato_slt,fact,zpo,1)
                     CALL calc_avg_xrms(avg_ca,avg_he,fstep,kavg_xrms
     &                    ,xpt0,ypt0,zpt0,xpo,ypo,zpo,wca,whe,wbc,protl
     &                    ,nato_slt)
                  END IF
               END IF
               
*=======================================================================
*---- Compute average RMS ----------------------------------------------
*=======================================================================
               
               IF(avg_rms) THEN
                  IF(avg_ca) CALL calc_avg2_str(protl,wca,xpt0,ypt0,zpt0
     &                 ,xp_avg2,yp_avg2,zp_avg2,xp0,yp0,zp0,nato_slt
     &                 ,iter_avg2)
                  IF(avg_he) CALL calc_avg2_str(protl,whe,xpt0,ypt0,zpt0
     &                 ,xp_avg2,yp_avg2,zp_avg2,xp0,yp0,zp0,nato_slt
     &                 ,iter_avg2)
               END IF
               IF(nrms .NE. 0) THEN
                  IF(MOD(nstep,nrms) .EQ.0) THEN
                     REWIND krms
                     CALL write_rms(krms,xp_avg,yp_avg,zp_avg,xp_avg2
     &                    ,yp_avg2,zp_avg2,iter_avg,fstep,ngrp,grppt
     &                    ,nres(1,1))
                  END IF 
               END IF
               
*=========== Compute dipole ============================================
               
               IF(ndipole.gt.0) THEN
                  IF(MOD(nstep,ndipole).EQ.0)THEN
                     CALL comp_dip(co,xpga,ypga,zpga,xpa,ypa,zpa,chrge
     &                    ,dips,ntap,ngrp,grppt)
                     WRITE(kdipole,106) fstep, (dips(j),j=1,3)
106                  FORMAT(' Dip. ',f12.3,3e15.5)
                  END IF
               END IF
               
               
*=========== Field and cofactors =======================================

               IF(ncofactor .NE. 0 .and.nvi.gt.0) THEN
                 IF(MOD(nstep,nvi) .EQ.0) THEN
                   aux = nstep*dtvi
                   CALL field_chromo(xpa,ypa,zpa,vi,cut_field,aux)
                 ENDIF
               ENDIF

c--            patchwork for computing torsional angle distribution 
c--            INACTIVE
c              call do_tors(xp0,yp0,zp0,ltor,ltors,betb) 
               
*=========== Plot .pdb file ============================================
               
               IF(nascii .NE. 0) THEN
                  IF(MOD(nstep,nascii) .EQ. 0) THEN
                     IF(ascii_nocell) THEN
                        CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xpo
     &                       ,ypo,zpo)
                     ELSE
                        CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot
     &                       ,protl)
                        CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo
     &                       ,ypo,zpo)
                     END IF
                     CALL plotc(fstep,beta,co,xpo,ypo,zpo,ntap,nres,m1
     &                    ,prsymb,idum)
                  END IF
               END IF
               
*=========== plot a fragment ===========================================
               
               IF(nplot_fragm .GT. 0) THEN
                  IF(MOD(nstep,nplot_fragm).EQ.0 ) THEN
                     write(kplot_fragm,*) ntot_fragm
                     write (kplot_fragm,*) 'fragments = ',nfragm 
                     do i=1,nfragm
                        fragm_1=fragm(1,i)
                        fragm_2=fragm(2,i)
                        CALL mts_plot_fragm(fragm_1,fragm_2,beta,xp0,yp0
     &                       ,zp0,ntap)
                     end do
                  END IF
               END IF
               
*========== Compute Hydrogen bonds =====================================
               
               IF(hbonds_tot .OR. hbonds_res) THEN
                  CALL write_hbonds(hbonds_tot,hbonds_res,fstep,khbonds
     &                 ,ss_index,nbun,nres,m1,co,xpa,ypa,zpa,lacc,ldon
     &                 ,llacc,lldon,rcut_hb,acut_hb,a2cut_hb)
               END IF
               
*========== Compute G of R's of solute and solvent =====================
               
               IF(gofr) THEN
                  IF(MOD(nstep,gofr_ncomp) .EQ. 0) THEN
                     l2=maxint
                     CALL calc_gofr(nato_slt,nato_slv,type_slv,ss_index
     &                    ,atomp,gofr_neighbor,gofr_intra,co,xpa,ypa,zpa
     &                    ,xpga,ypga,zpga,wca,whe,delrg,ntap,ngrp,grppt
     &                    ,l2,krdf,ngrdon,gofr_cut,nnlpp0)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
               END IF
               
               
*========== Print G of R's of solute and solvent =======================
               
               IF(gofr) THEN
                  IF(MOD(nstep,gofr_nprint) .EQ. 0)THEN
c                  vol_gofr=sum_volume/DFLOAT(nstep)
                     vol_gofr=volume
                     IF(slt_exist) THEN
                        offset=0
                        CALL write_gofrp(.NOT.gofr_avg,fstep,krdf,maxint
     &                       ,1,wca,whe,nato_slt,delrg,gofr_cut,ngrdon
     &                       ,iret,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                     IF(slv_exist .AND. slt_exist) THEN
                        CALL write_gofrw(.FALSE.,fstep,krdf,maxint,g1
     &                       ,delrg,gofr_cut,itype_slv,betab_slv
     &                       ,vol_gofr,ntype_slv,nmol,ngrdon,3,iret
     &                       ,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     ELSE IF(slv_exist .AND. (.NOT. slt_exist)) THEN
                        CALL write_gofrw(.TRUE.,fstep,krdf,maxint,g1
     &                       ,delrg,gofr_cut,itype_slv,betab_slv
     &                       ,vol_gofr,ntype_slv,nmol,ngrdon,3,iret
     &                       ,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                  END IF
               END IF
               
*=======================================================================
*---- Initialize G of Rs when required ---------------------------------
*=======================================================================
               
               IF(gofr .AND. gofr_avg) THEN
                  IF(MOD(nstep,gofr_navg) .EQ. 0) THEN
                     CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
                  END IF
               END IF
               
*=======================================================================
*----- Write topology --------------------------------------------------
*=======================================================================
               
               IF(prttopl) THEN
                  IF(MOD(nstep,ntop_print) .EQ. 0) THEN
                     IF(top_bonds(1) .GT. 0) CALL write_bonds(ktopol
     &                    ,fstep,top_bonds,lbnd,lbond,xp0,yp0,zp0)
                     IF(top_bendings(1) .GT. 0) CALL write_bends(ktopol
     &                    ,fstep,top_bendings,lbndg,lbend,xp0,yp0,zp0)
                     IF(top_ptors(1) .GT. 0) CALL write_tors('P',ktopol
     &                    ,fstep,top_ptors,ltor,ltors,xp0,yp0,zp0)
                     IF(top_itors(1) .GT. 0) CALL write_tors('I',ktopol
     &                    ,fstep,top_itors,litr,litor,xp0,yp0,zp0)
                  END IF
               END IF
               
*=======================================================================
*----- Write distance between fragments --------------------------------
*=======================================================================
               
               IF(fragm_dist) THEN
                  IF(MOD(nstep,nfragm_dist) .EQ. 0) THEN
                     CALL write_fragm_dist(fstep,kfragm_dist,fragm
     &                    ,nfragm,xpa,ypa,zpa,co)
                  END IF
               END IF
               
*=======================================================================
*----- Compute and write gyration ratio --------------------------------
*=======================================================================
               
               IF(wrtgyr) THEN
                  CALL write_gyr(kgyr,fstep,protl,nprot,ss_index,xp0,yp0
     &                 ,zp0)
               END IF
               
            ELSE
               nstep=nstep-1
            END IF
         END DO
      END IF
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp

*=======================================================================
*     Write timing
*=======================================================================

      WRITE(kprint,*)
      WRITE(kprint,60030)
      WRITE(kprint,17000) gcpu
      IF(time_corr) THEN
         WRITE(kprint,60300) gcpu/DFLOAT(ntap)
      ELSE
         WRITE(kprint,60200) gcpu/DFLOAT(nnstep)
      END IF
      WRITE(kprint,60030)

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(////20('*'),'  M o l e c u l a r   T o p o l o g y  ',
     x     20('*')/////)
1100  FORMAT('*',13(' '),
     x     ' I n p u t   O p e r a t i o n s   C o m p l e t e d ',
     x     12(' '),'*'/'*',78(' '),'*')
2000  FORMAT(//72('*')/'*',70(' '),'*'/
     &     '*     T R A J E C T O R Y   files contains  ',i6,
     &     '   configurations    *'/'*',70(' '),'*'/72('*')/)
2500  FORMAT(//72('*')/'*',70(' '),'*'/
     &     '*       C O M P U T I N G   velocity autocorrelation',
     &     ' function          *'/
     &     '*',70(' '),'*'/72('*')//)
2501  FORMAT(//72('*')/'*',70(' '),'*'/
     &     '*             C O M P U T I N G   atomic rms ',
     &     'displacements             *'/
     &     '*',70(' '),'*'/72('*')//)
3000  FORMAT(
     &     '                    ================================'/
     &     '                    =                              ='/
     &     '             -----> =  Processing atoms no. ',i5,'  ='/
     &     '                    =           through no. ',i5,'  ='/
     &     '                    =                              ='/
     &     '                    ================================'//)
3500  FORMAT(//'Writing v.a.c.f to file  --------> '//)
4500  FORMAT(//'Writing rms displacement to file  --------> '//)
5000  FORMAT(
     &     12x,'================================================='/
     &     12x,'=                                               ='/
     &     12x,'=                                               ='/
     &     12x,'=        SCANNING trajectory file ....          ='/
     &     12x,'=                                               ='/
     &     12x,'=                                               ='/
     &     12x,'================================================='/
     &     )
60030 FORMAT(/10x,'==========================================='/)
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
17000 FORMAT(  10x,' Total CPU time for data processing  = ',f10.3)
60200 FORMAT(  10x,' Averaged CPU time per step          = ',3x,f7.3)
60300 FORMAT(  10x,' Averaged CPU time per atomic corr.  = ',3x,f7.3)
90000 FORMAT(' Tstep ', f11.2,' fs successfully processed.',
     &' Program continues... ')
50000 FORMAT(10x,' Tstep =',f9.1)

80000 FORMAT('REMARK   Rigid body fit on CA atoms')
80100 FORMAT('REMARK   Rigid body fit on heavy atoms')
10977 FORMAT(//'*******WARNING: NO COFACTOR ATOMS SELECTED '/ 
     &     ' NCOFACTOR IS SET TO ZERO AND NO FIELD I COMPUTED'//)

      RETURN
      END



      SUBROUTINE do_tors(xp0,yp0,zp0,ltor,ltors,betb,fstep,dd)  

      IMPLICIT NONE

      REAL*8  xp0(*),yp0(*),zp0(*),fstep
      INTEGER ltors,ltor(4,*),it,l1,l2,l3,l4,i
      CHARACTER*7 betb(*)
      CHARACTER*1 b1,b2,b3,b4
      character*4 b44 

      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,xr4,yr4,zr4,gg1,dd
     &     ,dihed,gg2,dihed2
      real*8 x21,y21,z21,x32,y32,z32,x43,y43,z43,gx1,gy1,gz1,gx2,gy2,gz2

      DO i=1,ltors
          l1=ltor(1,i)
          l2=ltor(2,i)
          l3=ltor(3,i)
          l4=ltor(4,i)
          b1=betb(l1)
          b2=betb(l2)
          b3=betb(l3)
          b4=betb(l4)
          b44=b1//b2//b3//b4
          if(b44.ne.'cccc') go to 1000
          xr1=xp0(l1)
          yr1=yp0(l1)
          zr1=zp0(l1)
          xr2=xp0(l2)
          yr2=yp0(l2)
          zr2=zp0(l2)
          xr3=xp0(l3)
          yr3=yp0(l3)
          zr3=zp0(l3)
          xr4=xp0(l4)
          yr4=yp0(l4)
          zr4=zp0(l4)
          x21=xr2-xr1
          y21=yr2-yr1
          z21=zr2-zr1
          x32=xr3-xr2
          y32=yr3-yr2
          z32=zr3-zr2
          x43=xr4-xr3
          y43=yr4-yr3
          z43=zr4-zr3
          gx1 = y21*z32-z21*y32
          gy1 = z21*x32-x21*z32
          gz1 = x21*y32-y21*x32
          gx2 = y43*z32-z43*y32
          gy2 = z43*x32-x43*z32
          gz2 = x43*y32-y43*x32
          gg1 = dsqrt(gx1*gx1+gy1*gy1+gz1*gz1)
          gg2 = dsqrt(gx2*gx2+gy2*gy2+gz2*gz2)
          xr1 = (gy1*gz2-gz1*gy2)/(gg1*gg2)
          yr1 = (gz1*gx2-gx1*gz2)/(gg1*gg2)
          zr1 = (gx1*gy2-gy1*gx2)/(gg1*gg2)
          dihed = 180.d0*(1.d0 - dacos((gx1*gx2+gy1*gy2+gz1*gz2)/(gg1
     &         *gg2))/dacos(-1.d0))
          write(6,6554) b44,l1,l2,l3,l4,dihed,xr1,yr1,zr1
6554      FORMAT(a4,1x,4i4,3x,4f10.4) 
1000      CONTINUE
        END DO

        RETURN 
        end 
