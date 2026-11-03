      SUBROUTINE read_potential(fupdte,fabmd,err_args,err_unr,err_end)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Piero Procacci                                 *
*              Dip. Chimica,  Firenze                                  *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      REAL*8  fupdte,fabmd
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,nsevere,nwarning,j,kaux,i,ibeg,iend
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40),filepath
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      CHARACTER*10 chari
      LOGICAL  exist,ok1,ok2,ok3
      DATA sep/' ',','/comm/'#','!'/
      REAL*8 aux

*----------------------- EXECUTABLE STATEMENTS ------------------------*


c=======================================================================
c     Environment parser starts here 
c=======================================================================

      ok1=.false.
      ok2=.false.
      ok3=.false.
      j=0
      nat_added=0
      nat_removed=0
      nbonds_added = 0
      nbends_added = 0
      nitors_added = 0
      nbonds_steer_added = 0
      nbends_steer_added = 0
      nitors_steer_added = 0
      nsevere = 0 
      nwarning = 0 
      line(79:80)='  '
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

c==== Command GROUP_CUTOFF=============================================

      IF(strngs(1).EQ. 'GROUP_CUTOFF' ) THEN
         if(nword.eq.2.or.nword.eq.4) THEN 
            IF(strngs(2).EQ. 'ON' ) THEN
               grpcut=.TRUE.
               IF(nword .EQ. 4) THEN
                  IF(strngs(3) .EQ. 'CUTSPH') THEN
                     CALL fndfmt(2,strngs(4),fmt)
                     READ(strngs(4),fmt,err=20) rspon
                  END IF
               ELSE
                  nwarning = nwarning  + 1
                  errmsg='Unrecognized or missing keyword '// 
     &                 '; RSPON set to 1.0' 
                  CALL xerror(errmsg,80,1,11)
                  rspon=1.0D0
               END IF
            ELSE IF(strngs(2).EQ. 'OFF' ) THEN
               grpcut=.FALSE.
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3)//strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command EWALD====================================================

      ELSE IF(strngs(1).EQ. 'EWALD' ) THEN
         IF(nword .EQ. 1) THEN
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         ELSE
            IF(strngs(2).EQ. 'ON' ) THEN
               clewld=.TRUE.
               IF(nword .EQ. 4) THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) alphal
                  CALL fndfmt(2,strngs(4),fmt)
                  READ(strngs(4),fmt,err=20) rkcut
               ELSE IF(nword .EQ. 3) THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) alphal
               ELSE
                  nsevere = nsevere + 1
                  errmsg=err_args(1) // '1 after keyword "on"'
                  CALL xerror(errmsg,80,1,30)
               END IF
               
            ELSE IF(strngs(2).EQ. 'PME' ) THEN
               clewld=.true.
               pme=.true.
               IF(nword .EQ. 7) THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) alphal
                  CALL fndfmt(1,strngs(4),fmt)
                  READ(strngs(4),fmt,err=20) nfft1
                  CALL fndfmt(1,strngs(5),fmt)
                  READ(strngs(5),fmt,err=20) nfft2
                  CALL fndfmt(1,strngs(6),fmt)
                  READ(strngs(6),fmt,err=20) nfft3
                  CALL fndfmt(1,strngs(7),fmt)
                  READ(strngs(7),fmt,err=20) pme_order
               ELSE
                  nsevere = nsevere + 1
                  errmsg=err_args(1) // '4 after keyword "PME"' 
                  CALL xerror(errmsg,80,1,30)
               END IF
            ELSE IF(strngs(2) .EQ. 'REMOVE_MOMENTUM') THEN
               remove_momentum=.TRUE.
   
            ELSE IF(strngs(2).EQ. 'OFF' ) THEN
               clewld=.FALSE.
               
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3) // strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         END IF

c==== Command ERFC_SPLINE==============================================

      ELSE IF(strngs(1).EQ. 'ERFC_SPLINE' ) THEN
         erfc_spline=.TRUE.
         erfc_spline_corr=.false.
         IF(nword .NE. 1) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) erfc_bin
            IF(nword.NE.2) THEN 
              if(strngs(3).eq.'corrected') erfc_spline_corr=.true.
              if(nword.NE.3) THEN 
                CALL fndfmt(2,strngs(4),fmt)
                READ(strngs(4),fmt,err=20) rcut_corr
              END IF
            END IF
         END IF
         
c==== Command ERF_CORR=================================================

      ELSE IF(strngs(1).EQ. 'ERF_CORR' ) THEN
         erf_corr=.TRUE.
         IF(nword .EQ. 4) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nbinew
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rlew
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ruew
            delew = (ruew-rlew)/DFLOAT(nbinew) 
          ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) // '3 after keyword "ERF_CORR"' 
            CALL xerror(errmsg,80,1,30)
         ENDIF

c==== Command CUTOFF===================================================

      ELSE IF(strngs(1).EQ. 'CUTOFF' ) THEN
         IF(nword.ne.1) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) rspoff
         END IF
         
c==== Command LINKED CELL =============================================

      ELSE IF(strngs(1).EQ. 'LINKED_CELL' ) THEN
         LINKED_CELL = .TRUE.            
         IF(nword.EQ.4) THEN 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ncx
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ncy
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ncz
         ELSE IF(nword .EQ. 5) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ncx
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ncy
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ncz
            CALL fndfmt(1,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) nupdte_index
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'3'
            CALL xerror(errmsg,80,1,30)
         ENDIF

c==== Command VERLET_LIST==============================================

      ELSE IF(strngs(1).EQ. 'VERLET_LIST' ) THEN
         LINKED_CELL = .FALSE.            

c==== Command AGBNP ===================================================

      ELSE IF(strngs(1).EQ. 'AGBNP' ) THEN
         agbnp = .true.            

c==== Command UPDATE===================================================

      ELSE IF(strngs(1).EQ. 'UPDATE' ) THEN
         IF(nword .EQ. 3) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fupdte
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rspcut
         ELSE IF(nword .EQ. 5) THEN
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) hrcut
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) hacut
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'4'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  HBOND =================================================

      ELSEIF(strngs(1).EQ. 'H-BOND' ) THEN
         IF(strngs(2).EQ. 'ON' ) THEN
            hydbnd=.TRUE.
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) hrson
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) hrsoff
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) hanon
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) hanoff
            CALL fndfmt(1,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) nhskip
         ELSE IF(strngs(2).EQ. 'OFF' ) THEN
            hydbnd=.FALSE.
         ELSE
            errmsg=err_unr(3) // strngs(2) 
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  I-TORSION==============================================

      ELSE IF(strngs(1).EQ. 'I-TORSION' ) THEN
         IF( nword .NE. 2) THEN
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            IF(strngs(2) .EQ. 'HARMONIC' ) THEN
               itor_ptype=1
            ELSE IF(strngs(2) .EQ. 'COSINE' ) THEN
               itor_ptype=2
            ELSE
               errmsg=err_unr(3) // strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF
         
c==== Command  H-MASS=================================================

      ELSE IF(strngs(1).EQ. 'H-MASS' ) THEN
         hmass=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) hdmass
         
c==== Command  FFIELD =================================================

      ELSE IF(strngs(1).EQ. 'FFIELD' ) THEN
        call up_low(strngs(2),10)
        if(strngs(2)(1:5).eq."amber")  THEN 
          ffield_type="amber"
        else if(strngs(2)(1:4).eq."opls")  THEN 
          ffield_type="opls"
        end if

c==== Command  QQ-FUDGE===============================================
         
      ELSE IF(strngs(1) .EQ. 'QQ-FUDGE') THEN
         if(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fudge
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  LJ-FUDGE===============================================

      ELSE IF(strngs(1) .EQ. 'LJ-FUDGE') THEN
         if(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) lj_fudge
            lj_fudgeb=lj_fudge
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  BENDING================================================

      ELSE IF(strngs(1) .EQ. 'BENDING') THEN 
         if(nword.eq.2) THEN 
            IF(strngs(2) .EQ. 'OFF') THEN
               bending=.FALSE.
            ELSE IF(strngs(2) .EQ. 'ON') THEN
               bending=.TRUE.
            ELSE 
               errmsg=err_unr(3)//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
            
c==== Command  AUTO_DIHEDRAL==========================================

      ELSE IF(strngs(1).EQ. 'AUTO_DIHEDRAL' ) THEN
         adihed = .TRUE.
         
c==== Command  SELECT_DIHEDRAL========================================

      ELSE IF(strngs(1).EQ. 'SELECT_DIHEDRAL' ) THEN
         adihed = .FALSE.
         
c==== Command  SELECT_STRETCHING======================================

      ELSE IF(strngs(1).EQ. 'STRETCHING' ) THEN
         stretch = .TRUE.
         IF(nword.gt.1) THEN 
            IF(strngs(2) .EQ. 'HEAVY') THEN
               stretch_heavy= .TRUE.
            ELSE 
               errmsg=err_unr(3)//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF

c==== Command  COMMAND NONBONDED_OFF==================================

      ELSE IF(strngs(1).EQ. 'NONBONDED_OFF' ) THEN
         nonbnd=.FALSE.
         
c==== Command  COMMAND NONBONDED_OFF==================================

      ELSE IF(strngs(1).EQ. 'CHECK_COORD_OFF' ) THEN
         check_coord=.false.
         
c==== Command  CONSTRAINT  ===========================================

      ELSE IF(strngs(1).EQ. 'CONSTRAINT') THEN
         IF(nword .GT. 1) THEN 
            IF(strngs(2) .EQ. 'SHAKE') THEN
               mim_lim=0
            ELSE IF(strngs(2) .EQ. 'MIM') THEN
               IF(nword .EQ. 3) THEN
                  CALL fndfmt(1,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) mim_lim
               ELSE
                  mim_lim=20
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3)//strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  KEEP_BONDS =============================================

      ELSE IF(strngs(1).EQ. 'KEEP_BONDS' ) THEN
         adjust_cnstr = .FALSE.

c==== Command  DEFINE_ALCHEMICAL_ATOMS ========================================

      ELSE IF(strngs(1).EQ. 'DEFINE_ALCHEMICAL_ATOM' ) THEN
          if(path_steer) THEN
           errmsg = 'cannot define alchemical atoms *after* STEER_PATH'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
         if(nat_added.ge.1000.OR.nat_removed.GT.1000) THEN 
           errmsg =
     &    ' DEFINE ALCHEMICAL ATOM: too many added or removed atoms :'
     &          // ' Max is 1000'
           CALL xerror(errmsg,80,1,30)
           STOP
         END IF
         if(nword.eq.4) THEN  
           if(strngs(4).eq."on") THEN
             CALL fndfmt(1,strngs(2),fmt)
             READ(strngs(2),fmt,err=20)  ibeg 
             CALL fndfmt(1,strngs(3),fmt)
             READ(strngs(3),fmt,err=20)  iend 
             do i=ibeg,iend
               nat_added = nat_added+1
               atom_added(nat_added)=i
             end do
           else if (strngs(4).eq."off") THEN
             CALL fndfmt(1,strngs(2),fmt)
             READ(strngs(2),fmt,err=20)  ibeg 
             CALL fndfmt(1,strngs(3),fmt)
             READ(strngs(3),fmt,err=20)  iend 
             do i=ibeg,iend
               nat_removed = nat_removed+1
               atom_removed(nat_removed)=i
             end do
           else 
             nsevere = nsevere + 1
             errmsg= " Use 'on/off' to define alchemical switching" 
             CALL xerror(errmsg,80,1,30)
           end if
         ELSE
           errmsg=err_args(1) // '4'
           CALL xerror(errmsg,80,1,30)
         END  IF


c== new alchemical command
!       ALCHEMY 
!         lig_on ibeg1 iend1 (name_on)
!         lig_on ibeg2 iend2 (name_on)
!         ..                       
!         lig_off ibeg1 iend1 (name_off) 
!         lig_off ibeg2 iend2 (name_off)
!         ..                       
!           (read lig off protocols from name_off) 
!           (read lig on protocols from name_on) 
!         res_on  ibeg1 iend1 Cbeta (res_name_on) 
!         res_on  ibeg2 iend2 Cbeta (res_name_on) 
!         .. 
!         res_off  ibeg1 iend1 Calpha (res_name_on) 
!         res_off  ibeg2 iend2 Calpha (res_name_on) 
!           (read res on protocols from res_name_on)
!           (read res off  protocols from res_name_ff)
!           (attach CB to CA (and add all the rest of topology)           
!           NB PDB carries on top on species.
!           very careful for off species (ibeg, iend correct)
!           do a careful printout in output to show that
!           you did the right thing.  
!        END

c==== Command  ADD_BONDS ========================================

      ELSE IF(strngs(1).EQ. 'ADD_STR_BONDS' ) THEN
         nbonds_added = nbonds_added+1
         strbonds(nbonds_added)=0
         if(path_steer) THEN
           errmsg = 'cannot define extra steer coord *after* STEER_PATH'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
         if(nbonds_added.ge.500) THEN 
           errmsg =
     &          ' toomany ADD_STR_BONDS directives:'
     &          // ' Max is 500'
           CALL xerror(errmsg,80,1,30)
           STOP
         END IF
         if(nword.eq.5.or.nword.eq.6) THEN  
           CALL fndfmt(1,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) atom_b1(nbonds_added)
           CALL fndfmt(1,strngs(3),fmt)
           READ(strngs(3),fmt,err=20) atom_b2(nbonds_added)
           if(atom_b1(nbonds_added).eq.atom_b2(nbonds_added)) THEN 
             errmsg = ' atom 1 must be not equal to atom 2 ' 
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
           CALL fndfmt(2,strngs(4),fmt)
           READ(strngs(4),fmt,err=20) force_cost(nbonds_added)
           CALL fndfmt(2,strngs(5),fmt)
           READ(strngs(5),fmt,err=20) eqdist(nbonds_added)
           IF(nword.eq.6) THEN 
             CALL fndfmt(2,strngs(6),fmt)
             READ(strngs(6),fmt,err=20) eqdist1(nbonds_added)
             aux = abs(eqdist1(nbonds_added)-eqdist(nbonds_added))
             if(aux.gt.0.001D0.and.eqdist1(nbonds_added).gt.0) THEN  
               strbonds(nbonds_added) = 1
             else
               strbonds(nbonds_added) = 1 ! K switching on : not active yet
             end if
             nbonds_steer_added=nbonds_steer_added+1
           END IF
         ELSE
           nsevere = nsevere + 1
           errmsg=err_args(1) // '5 or 6'
           CALL xerror(errmsg,80,1,30)
         END  IF

c==== Command  ADD_BENDINGS ===========================================

      ELSE IF(strngs(1).EQ. 'ADD_STR_BENDS' ) THEN
         nbends_added = nbends_added+1
         strbends(nbends_added)=0
         if(path_steer) THEN
           errmsg = 'cannot define extra steer coord *after* STEER_PATH'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
         if(nbends_added.ge.500) THEN 
           errmsg =
     &          ' toomany ADD_STR_BENDS directives:'
     &          // ' Max is 500'
           CALL xerror(errmsg,80,1,30)
           STOP
         END IF
         if(nword.eq.6.or.nword.eq.7) THEN  
           CALL fndfmt(1,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) atom_be1(nbends_added)
           CALL fndfmt(1,strngs(3),fmt)
           READ(strngs(3),fmt,err=20) atom_be2(nbends_added)
           CALL fndfmt(1,strngs(4),fmt)
           READ(strngs(4),fmt,err=20) atom_be3(nbends_added)

           if(atom_be1(nbends_added).eq.atom_be2(nbends_added).or. 
     &        atom_be2(nbends_added).eq.atom_be3(nbends_added).or. 
     &        atom_be1(nbends_added).eq.atom_be3(nbends_added))
     &     THEN 
             errmsg = ' atom labels must be different '
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
           CALL fndfmt(2,strngs(5),fmt)
           READ(strngs(5),fmt,err=20) force_ang(nbends_added)
           CALL fndfmt(2,strngs(6),fmt)
           READ(strngs(6),fmt,err=20) eqang(nbends_added)
           IF(nword.eq.7)THEN 
             CALL fndfmt(2,strngs(7),fmt)
             READ(strngs(7),fmt,err=20) eqang1(nbends_added)
             aux = abs(eqang1(nbends_added)-eqang(nbends_added))
             if(aux.gt.0.001D0.and.eqang1(nbends_added).lt.180.d0) THEN  
               strbends(nbends_added) = 1
             else
               strbends(nbends_added) = 1 ! K switching: not active yet 
             end if
             nbends_steer_added=nbends_steer_added+1
           END IF
         ELSE
           nsevere = nsevere + 1
           errmsg=err_args(1) //'6 or 7'
           CALL xerror(errmsg,80,1,30)
         END  IF

c==== Command  ADD_TORSION ===========================================

      ELSE IF(strngs(1).EQ. 'ADD_STR_TORS' ) THEN
         nitors_added = nitors_added+1
         strtors(nitors_added)=0
         if(path_steer) THEN
           errmsg = 'cannot define extra steer coord *after* STEER_PATH'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
           
         if(nitors_added.ge.500) THEN 
           errmsg =
     &          ' toomany ADD_STR_TORS directives:'
     &          // ' Max is 500'
           CALL xerror(errmsg,80,1,30)
           STOP
         END IF
         if(nword.eq.7.or.nword.eq.8) THEN  
           CALL fndfmt(1,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) atom_it1(nitors_added)
           CALL fndfmt(1,strngs(3),fmt)
           READ(strngs(3),fmt,err=20) atom_it2(nitors_added)
           CALL fndfmt(1,strngs(4),fmt)
           READ(strngs(4),fmt,err=20) atom_it3(nitors_added)
           CALL fndfmt(1,strngs(5),fmt)
           READ(strngs(5),fmt,err=20) atom_it4(nitors_added)

           if(atom_it1(nitors_added).eq.atom_it2(nitors_added).or. 
     &        atom_it1(nitors_added).eq.atom_it3(nitors_added).or. 
     &        atom_it1(nitors_added).eq.atom_it4(nitors_added).or. 
     &        atom_it2(nitors_added).eq.atom_it3(nitors_added).or. 
     &        atom_it2(nitors_added).eq.atom_it4(nitors_added).or. 
     &        atom_it3(nitors_added).eq.atom_it4(nitors_added))
     &     THEN 
             errmsg = ' atom labels must be different '
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
           CALL fndfmt(2,strngs(6),fmt)
           READ(strngs(6),fmt,err=20) force_died(nitors_added)
           CALL fndfmt(2,strngs(7),fmt)
           READ(strngs(7),fmt,err=20) eqdied(nitors_added)
           if(eqdied(nitors_added).lt.-180.d0) THEN 
             errmsg = ' init diehdral angle cannot be  lt -180 deg. '
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
           END IF
           if(nword.eq.8) THEN 
             CALL fndfmt(2,strngs(8),fmt)
             READ(strngs(8),fmt,err=20) eqdied1(nitors_added)
             aux = abs(eqdied1(nitors_added)-eqdied(nitors_added))
             if(aux.gt.0.001D0.and.eqdied1(nitors_added).lt.360.d0) THEN  
               strtors(nitors_added) = 1
             else
               strtors(nitors_added) = 1 ! K switching: on not active yet 
             end if
             nitors_steer_added=nitors_steer_added+1
           END IF
         ELSE
           nsevere = nsevere + 1
           errmsg=err_args(1) //'8'
           CALL xerror(errmsg,80,1,30)
         END  IF

c==== Command  STEER_PATH  ==========================================

       ELSE IF(strngs(1).EQ. 'STEER_PATH' ) THEN
         path_steer=.true.
         IF(  nbonds_steer_added.EQ.0.and.
     &        nbends_steer_added.EQ.0.and.
     &        nitors_steer_added.EQ.0.and.
     &        nat_added.EQ.0.and.nat_removed.EQ.0) THEN
           errmsg='No steering coordinates or alchemical atoms defined.'
     &         //' Define them *before* STEER_PATH'  
           CALL xerror(errmsg,80,1,30)
           nsevere=nsevere+1
         END IF
         if(nbonds_steer_added.gt.500) THEN
           errmsg='More than 500 bonds steer defined' 
           CALL xerror(errmsg,80,1,30)
           nsevere=nsevere+1
         end if  
         if(nbends_steer_added.gt.500) THEN
           errmsg='More than 500 bends steer defined' 
           CALL xerror(errmsg,80,1,30)
           nsevere=nsevere+1
         end if  
         if(nitors_steer_added.gt.500) THEN
           errmsg='More than 500 tors steer defined' 
           CALL xerror(errmsg,80,1,30)
           nsevere=nsevere+1
         end if  
         IF(strngs(2) .EQ. 'OPEN') THEN
           CALL uscrpl(strngs(4),80)
           filepath=strngs(3)
           INQUIRE(FILE=filepath,EXIST=exist)
           IF(.NOT. exist) THEN
             errmsg=
     &         'Path File '//TRIM(filepath)//' does not exist. Abort.'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
           ELSE
             CALL openf(kaux,strngs(3),'FORMATTED','OLD',0)
             do i=1,100
               read(kaux,*,end=1008,ERR=1007)  timerec(i),
     &              (pathbo(i,j),activebo(i,j),j=1,nbonds_steer_added),
     &              (pathbe(i,j),activebe(i,j),j=1,nbends_steer_added),
     &              (pathto(i,j),activeto(i,j),j=1,nitors_steer_added)
               timerec(i)=timerec(i)*1000.d0
             end do
             write(chari,'(f10.2)') timerec(100) 
             errmsg='Path file too long; Last timeslice read in '//chari
!            reset steered X coord to type 1 if 6th arg of ADD_STR_X  was
!            negative (K switch); no k switch in steer path. 
             do i=1,nbonds_added
               if(strbonds(i).EQ.2) strbonds(i)=1 
             end do
             do i=1,nbends_added
               if(strbends(i).EQ.2) strbends(i)=1 
             end do
             do i=1,nitors_added
               if(strtors(i).EQ.2) strtors(i)=1 
             end do

             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere +1 
1007         errmsg=
     &            'Something wrong in reading Path File'//TRIM(filepath)
             CALL xerror(errmsg,80,1,30)
             write(kprint,1009) nbonds_steer_added,nbends_steer_added
     &            ,nitors_steer_added
1009         format(' ===== On each line ORAC expects: time ',i4
     &            ,' bonds',i4,' bends and',i4,' tors'/
     &            ' ===== where time must be expressed in picoseconds')
             nsevere=nsevere+1
1008         ntimes=i-1
           END IF
c   Path steer for alchemical types
         ELSE IF(strngs(2) .EQ. 'ALCHEMY') THEN
           alchemy=.true.
           CALL uscrpl(strngs(4),80)
           filepath=strngs(3)
           INQUIRE(FILE=filepath,EXIST=exist)
           IF(.NOT. exist) THEN
             errmsg=
     &         'Alchemical File '//TRIM(filepath)// 
     &            ' does not exist. Abort.'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
           ELSE
             CALL openf(kaux,strngs(3),'FORMATTED','OLD',0)
             do i=1,100
c              read time protocol for added and removed atoms
               read(kaux,*,end=1018,ERR=1017)  timerec(i),ladded(i)
     &              ,laddedq(i),lremoved(i),lremovedq(i)
               timerec(i)=timerec(i)*1000.d0
             end do
1017         write(chari,'(f10.2)') timerec(i) 
             errmsg='Alchemical file wrong or too long; Last timeslice'
     &           //'read in'//chari
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere +1 
1018         ntimes=i-1
           END IF
         ELSE
           errmsg='OPEN keyword not found'
           CALL xerror(errmsg,80,1,30)
           nsevere = nsevere + 1
         END IF
        
!==== Command  ADD_STR_COM  ==========================================

      ELSE IF(strngs(1).EQ. 'ADD_STR_COM' ) THEN
c------- read the line
         addstrcom=.true.
 155     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 155
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         
c------  subcommand "ligand"---------------------------------------------- 
         IF(strngs(1).EQ.'ligand') THEN
            CALL fndfmt(1,strngs(2),fmt) 
            READ(strngs(2),fmt,err=20) ilig1
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ilig2 
            do i=ilig1,ilig2
              indlig(i)=.true. 
            end do
            ok1=.true.

c------  subcommand "target"---------------------------------------------- 
         ELSE IF(strngs(1).EQ.'target') THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20)  itar1
            CALL fndfmt(1,strngs(3),fmt) 
            READ(strngs(3),fmt,err=20)  itar2
            do i=itar1,itar2
              indtar(i)=.true. 
            end do
            ok2=.true.

c------  subcommand "force" 
         ELSE IF(strngs(1).EQ.'force') THEN
           if(nword.eq.3) THEN 
             CALL fndfmt(2,strngs(2),fmt)
             READ(strngs(2),fmt,err=20) kcom
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) r0com
             ok3=.true.
           else  if(nword.eq.4) THEN
             steer_com=.true.
             CALL fndfmt(2,strngs(2),fmt)
             READ(strngs(2),fmt,err=20) kcom
             CALL fndfmt(2,strngs(3),fmt)
             READ(strngs(3),fmt,err=20) r0com
             CALL fndfmt(2,strngs(4),fmt)
             READ(strngs(4),fmt,err=20) r1com
             ok3=.true.
           else 
             nsevere = nsevere + 1
             errmsg=err_args(1) //'2 or 3'
             CALL xerror(errmsg,80,1,30)
           end if
         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1).EQ. 'END' ) THEN
           if(.not.ok1) THEN 
             nsevere = nsevere + 1
             errmsg="ADD_STR_COM: ligand not provided"  
             CALL xerror(errmsg,80,1,30)
           endif
           if(.not.ok2) THEN 
             nsevere = nsevere + 1
             errmsg="ADD_STR_COM: target not provided"  
             CALL xerror(errmsg,80,1,30)
           endif
           if(.not.ok3) THEN 
             nsevere = nsevere + 1
             errmsg="ADD_STR_COM: k and r0 not provided"  
             CALL xerror(errmsg,80,1,30)
           endif
           GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF
         GOTO 155

c==== Command  ADJUST_BONDS ==========================================

      ELSE IF(strngs(1).EQ. 'ADJUST_BONDS' ) THEN
        adjust_cnstr = .TRUE.

c==== Command  INTRAMOL      ==========================================

      ELSE IF(strngs(1).EQ. 'INTRAMOL' ) THEN
         intramol = .TRUE.
         if(nword.EQ.2) THEN
           if(strngs(2).EQ."charge".OR.strngs(2).EQ."CHARGE") THEN
             ljintra=.true.
           ELSE
c---        could not fine SUBCOMMAND charge/CHARGE
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
          END IF
        ENDIF
c==== Command  SOFT-CORE  ================alchemy==========================

      ELSE IF(strngs(1).EQ. 'SOFT-CORE' ) THEN
         shifted = .FALSE.
         IF(nword.EQ.3) THEN
           CALL fndfmt(2,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) qsoft
           CALL fndfmt(2,strngs(3),fmt)
           READ(strngs(3),fmt,err=20) ljsoft
         ELSE
           nsevere = nsevere + 1
           errmsg=err_args(1) //'3'
           CALL xerror(errmsg,80,1,30)
         END  IF

c==== Command  SHIFT-POT  ===================alchemy=======================

      ELSE IF(strngs(1).EQ. 'SHIFT-POT' ) THEN
         shifted = .TRUE.
         IF(nword.EQ.3) THEN
           CALL fndfmt(2,strngs(2),fmt)
           READ(strngs(2),fmt,err=20) qshift
           CALL fndfmt(2,strngs(3),fmt)
           READ(strngs(3),fmt,err=20) ljshift
         ELSE
           nsevere = nsevere + 1
           errmsg=err_args(1) //'3'
           CALL xerror(errmsg,80,1,30)
         END  IF
c==== Command BLANK LINE===============================================

      ELSE IF(strngs(1).EQ. ' ') THEN

c==== Begininning of next ENVIRONMENT =================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
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
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_POTENTIAL'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,1)
      ENDIF    

c==============================================================================
c     Verification Part
c==============================================================================

c--   Check if pme_orded is gt 3 
      if(pme.and.(pme_order.lt.3)) THEN 
         errmsg='"Ewald": PME order must be at least 3'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF
      if(linked_cell) THEN 
         if(ncx.le.0.or.ncy.le.0.or.ncz.le.0) THEN 
            errmsg='Bad Linked Cell Neighbor list parameters'  
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
      END IF

      if(erf_corr.and.alchemy) THEN 
        errmsg="ERF_CORR used for ALCHEMY calculations" 
        CALL xerror(errmsg,80,1,1)
        errmsg="Be sure that interval includes all intra contacts" 
        CALL xerror(errmsg,80,1,1)
      ENDIF

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_POTENTIAL'
         CALL xerror(errmsg,80,1,2)
         STOP 
      ELSE IF(nsevere.gt.99) THEN
         errmsg='MORE THAN 99 ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      END IF

      RETURN

c==============================================================================
c     EOF found while reading ABMD subcommands
c==============================================================================

300   CONTINUE
      iret=1
      errmsg='EOF found while reading ABMD subcommands'
      CALL xerror(errmsg,80,1,2)
      STOP

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
