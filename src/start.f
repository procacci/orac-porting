      SUBROUTINE starta(mapnl,mapnl_slv,mapp,mapp_slv)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine is called once at the beginning of the run.      *
*     It calls other routines which read input data for the run and    *
*     perform the initialisation of many variables. There are 11       *
*     input arguments to be provided which define the physical         *
*     dimensions of the arrays that are passed.                        *
*                                                                      *
*                                                                      *
*   Time-stamp: <97/02/07 15:46:56 marchi>                             *
*                                                                      *
*     Written by Massimo Marchi CEA Saclay 1997                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*    STARTA externals: 	       	       	       	       	       	       *
*        brot change_tors chmass chnpr comp_concth		       *
*        covbd covchg covfod covnb covpr covuen			       *
*        dumptp dustar igmap join mappa min_pack		       *
*        pntres readco readtp read_input redmss rmpbond		       *
*        scale_charges search_clsth setuns set_ss_array spec14	       *
*        xerror							       *
*       							       
************************************************************************


*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER mapnl(*),mapnl_slv(*)
      LOGICAL   mapp(*),mapp_slv(*)

*=======================================================================
*     NORES   :  Maximum number of residue units.
*=======================================================================

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER l1,iret,ksitew,nmapnl,nmapnl_slv
      CHARACTER*80 errmsg
      REAL*8 hana,hanb
      LOGICAL near0,mask(m1),mask2(m1),slv_exist_old
      EXTERNAL near0
      REAL*8  useed
      REAL*8  dustar,dunib,ubits
      CHARACTER*9 mesg


*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*----- Read in all DATA needed at running time -------------------------
*=======================================================================

      CALL read_input(ksitew,iret,errmsg)

*----------------- If iret .EQ. 1 STOP !! ------------------------------
*                                 ====

          IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*----------------- If iret .EQ. 2 a WARNING is issued ------------------
*                                   =======

*
*
          IF(iret.EQ.2) CALL xerror(errmsg,80,1,1)


*=======================================================================
*----- Initialize seed in random number generator and do test ----------
*----- of the generator ------------------------------------------------
*=======================================================================

      useed=dustar(iseed)


*=======================================================================
*----- Setup working units of the program ------------------------------
*=======================================================================

      CALL setuns(time)

*=======================================================================
*--------- Compute the box to orthogonal frame matrix and its inverse --
*=======================================================================

      IF(cpress .AND. (nflag(1) .NE. 0)) THEN
         CALL readco(co,oc,volume)
      END IF
      IF(nflag(1).EQ.0) THEN 
        IF(.NOT. stoprun .AND. .NOT. analys) THEN
          CALL brot(read_co,aaxis,baxis,caxis,alf,bet,gam,icl,icm,icn
     &         ,co,oc,volume)
        END IF
      END IF
*+++++++++++++++ INITIALISE MOLECULE +++++++++++++++++++++++++++++++++++


*=======================================================================
*----------- Change improper torsion parametes and arrays --------------
*=======================================================================

      CALL change_tors(itor_ptype)

      IF(tpgfil) THEN

*=======================================================================
*----- Read in the topological file ------------------------------------
*=======================================================================
        write(kprint,1000)
1000    FORMAT(//5x,'<------ Reading Topology-Parameter File ------->'
     &       //)
         CALL readtp(ktpgprm_read,mapnl,nmapnl)
         IF(slv_create .OR. slv_add) THEN
            IF(slv_exist) THEN
               errmsg=
     &'Solvent was found in the TPGPRM file: Cannot create'
     &//' or add more.'
               CALL xerror(errmsg,80,1,2)
            END IF
            nmol=nmol_create
            slv_exist=.TRUE.
            mesg='SOLVENT:  '
            CALL join(nbun_slv,nato_slv,ngrp_slv,lbond_slv,lbend_slv
     &           ,ltors_slv,litor_slv,lphyd,lpnbd,nbone_slv
     &           ,int14p_slv,int13p_slv,llacc_slv,lldon_slv,mback_slv
     &           ,grppt_slv,int14_slv,int13_slv,nres_slv,mend_slv
     &           ,lacc_slv,ldon_slv,nhtype_slv,lbnd_slv,lbndg_slv
     &           ,ltor_slv,litr_slv,atres_slv,concta_slv,nbtype_slv
     &           ,mass_slv,beta_slv,betb_slv,alnbd,potbo_slv
     &           ,potbe_slv,potto_slv,potit_slv,chrge_slv,bsitp_slv
     &           ,asitp_slv,nbsitp_slv,nasitp_slv,nrigg_slv,prsymb
     &           ,slvatm,slv3,slv4,slv5,types,slv2,slv6,slv_group
     &           ,slv_cns1,slv_cns2,debug_bt,debug_pt,debug_it,debug_rs
     &           ,debug_ct,debug_st,adihed,mesg)

*=======================================================================
*--------- Create a list of interaction types for third neighbour ------
*--------- interactions ------------------------------------------------
*=======================================================================
            
            CALL spec14(int14_slv,int14p_slv,betb_slv,nbtype_slv
     &           ,alnbd,lpnbd,nato_slv,type14_slv)

*=======================================================================
*--------- Produce the non-bonded interactions map. Add ----------------
*--------- to MAPNL addresses for 1-2, 1-3 and 1-4 ---------------------
*=======================================================================

            CALL mappa(lbnd_slv,lbond_slv,lbndg_slv,lbend_slv,int14_slv
     &           ,int14p_slv,nato_slv,mapnl_slv,nmap_slv,mapp_slv)


*=======================================================================
*----- Write the topological file when required ------------------------
*=======================================================================

            CALL min_pack(nato_slv,mapnl_slv,nmapnl_slv)

         END IF
      ELSE

************************************************************************
******* DO SOLUTE FIRST ************************************************
************************************************************************

*=======================================================================
*----------- Assemble the solute and define its topology ---------------
*=======================================================================
      
         IF(slt_exist) THEN
            mesg='SOLUTE:  '
            CALL join(nbun,ntap,ngrp,lbond,lbend,ltors,litor,lphyd,lpnbd
     &           ,nbone,int14p,int13p,llacc,lldon,mback,grppt,int14
     &           ,int13,nres,mend,lacc,ldon,nhtype,lbnd,lbndg,ltor,litr
     &           ,atres,concta,nbtype,mass,beta,betb,alnbd,potbo,potbe
     &           ,potto,potit,chrge,bsitp,asitp,nbsitp,nasitp,nrigg
     &           ,prsymb,m1,m2,m3,m4,m6,m9,m10,m11,m13,m14,debug_bt
     &           ,debug_pt,debug_it,debug_rs,debug_ct,debug_st,adihed
     &           ,mesg)

*=======================================================================
*--------- Create a list of interaction types for third neighbour ------
*--------- interactions ------------------------------------------------
*=======================================================================

            CALL spec14(int14,int14p,betb,nbtype,alnbd,lpnbd,ntap,type14
     &           )

*=======================================================================
*--------- Produce the non-bonded interactions map. Add ----------------
*--------- to MAPNL addresses for 1-2, 1-3 and 1-4 ---------------------
*=======================================================================

            CALL mappa(lbnd,lbond,lbndg,lbend,int14,int14p,ntap,mapnl,
     &        m8,mapp)


*=======================================================================
*----- Write the topological file when required ------------------------
*=======================================================================

            CALL min_pack(ntap,mapnl,nmapnl)
         END IF

************************************************************************
******* DO THEN THE SOLVENT IF REQUIRED ********************************
************************************************************************

*=======================================================================
*----------- Assemble the solvent and define its topology --------------
*=======================================================================
         IF(slv_create .OR. slv_add) THEN
            slv_exist=.TRUE.
            nmol=nmol_create
            mesg='SOLVENT: '
            CALL join(nbun_slv,nato_slv,ngrp_slv,lbond_slv,lbend_slv
     &           ,ltors_slv,litor_slv,lphyd,lpnbd,nbone_slv
     &           ,int14p_slv,int13p_slv,llacc_slv,lldon_slv,mback_slv
     &           ,grppt_slv,int14_slv,int13_slv,nres_slv,mend_slv
     &           ,lacc_slv,ldon_slv,nhtype_slv,lbnd_slv,lbndg_slv
     &           ,ltor_slv,litr_slv,atres_slv,concta_slv,nbtype_slv
     &           ,mass_slv,beta_slv,betb_slv,alnbd,potbo_slv
     &           ,potbe_slv,potto_slv,potit_slv,chrge_slv,bsitp_slv
     &           ,asitp_slv,nbsitp_slv,nasitp_slv,nrigg_slv,prsymb
     &           ,slvatm,slv3,slv4,slv5,types,slv2,slv6,slv_group
     &           ,slv_cns1,slv_cns2,debug_bt,debug_pt,debug_it,debug_rs
     &           ,debug_ct,debug_st,adihed,mesg)
            
*=======================================================================
*--------- Create a list of interaction types for third neighbour ------
*--------- interactions ------------------------------------------------
*=======================================================================
            
            CALL spec14(int14_slv,int14p_slv,betb_slv,nbtype_slv
     &           ,alnbd,lpnbd,nato_slv,type14_slv)
         
*=======================================================================
*--------- Produce the non-bonded interactions map. Add ----------------
*--------- to MAPNL addresses for 1-2, 1-3 and 1-4 ---------------------
*=======================================================================
            
            CALL mappa(lbnd_slv,lbond_slv,lbndg_slv,lbend_slv,int14_slv
     &           ,int14p_slv,nato_slv,mapnl_slv,nmap_slv,mapp_slv)
         

*=======================================================================
*----- Write the topological file when required ------------------------
*=======================================================================

            CALL min_pack(nato_slv,mapnl_slv,nmapnl_slv)
         END IF
         IF(slt_exist .AND. .NOT. (slv_create .OR. slv_add)) THEN
            CALL set_ss_array(ss_point,m1+1,ss_index,nmol,nato_slv,ntap)
            IF(tpgwbn .AND. (.NOT. sgroup)) THEN
               CALL dumptp(ktpgprm_write,mapnl,nmapnl)
            END IF
         END IF
      END IF

*=======================================================================
*----- Change charges to program units ---------------------------------
*=======================================================================

      IF(nflag(1) .NE. 0 .AND. sgroup) THEN
         errmsg=
     &'In START: Symmetry operations cannot be carried out with '
     &            //' CONTROL .NE. 0. Abort.'
         CALL xerror(errmsg,80,1,2)
      END IF

      IF(.NOT. sgroup .AND. (.NOT. slv_create) .AND. (.NOT. slv_add))
     &     THEN
*=======================================================================
*----  Decide which stretching and constraint to keep ------------------
*=======================================================================

         CALL rmpbond(stretch,stretch_heavy,lbnd,lbond,lstrtch
     &        ,lstretch,lcnstr,lconstr,M9,betb,potbo,potbo_cnst)

*=======================================================================
*----  Compute a new connection table including the new constraints ----
*=======================================================================

         IF(stretch_heavy) THEN
            CALL comp_concth(ntap,beta,concta,concth,M1)
            L1=2*m1
            CALL search_clsth(clsthl,L1,nclsth,concth,m1,mask,mask2,ntap
     &           ,iret,errmsg)
            IF(iret.NE.0) RETURN
         END IF

*=======================================================================
*-------- Find out how many molecules form the solute ------------------
*=======================================================================

         L1=2*m1
         CALL chnpr(protl,L1,npm,nprot,concta,m1,mask,mask2,ntap,iret
     &        ,errmsg)
         IF(iret .NE. 0) THEN
            CALL xerror(errmsg,80,1,2)
         END IF

*=======================================================================
*----- Compute the total charge and scale it to zero when required -----
*=======================================================================

         CALL scale_charges(kprint,nprot_charges,prot_charges
     &        ,chrge,protl,nprot,scharge)

         CALL covchg(chrge,ntap)
             
*=======================================================================
*--------- Change the covalent interaction parameters ------------------
*---------------- program units ----------------------------------------
*=======================================================================

         CALL covbd(potbo,potbe,potto,potit,ptorj,lstretch,lbend
     &            ,ltors,litor,itor_ptype,m9,m2,m3,m4)

*=======================================================================
*--------- Change the non bonded interaction parameters to program -----
*--------- units and compute the macromolecule total mass --------------
*=======================================================================

         CALL covnb(pnbd1,pnbd2,pnbd3,pnbd4,sjorg,ejorg,lpnbd,iz
     &            ,phyd1,phyd2,lphyd,mass,wmtp,ntap,ecc6,ecc12,ecc146
     &            ,ecc1412,c6jorg,c12jorg,type_table,lj_fudge,m6)

*=======================================================================
*--------- Compute the intergroup interaction map ----------------------
*=======================================================================

         CALL igmap(ngrp,grppt,ingrpp,ingrp,m12,mapnl,errmsg,iret)
         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--------- Compute pointers to start an end of residue -----------------
*=======================================================================

         CALL pntres(nres(1,1),mres,atomg,grppt,nbun,resg,ngrp,ntap)
      END IF
          
*=======================================================================
*----- Change the mass of the hydrogens when required ------------------
*=======================================================================

      IF(hmass) CALL chmass(beta,mass,ntap,hdmass)

*=======================================================================
*--------- Compute reduced mass for atomic groups ----------------------
*=======================================================================

      CALL redmss(ngrp,grppt,mass,pmass)
          
*=======================================================================
*----- Change the ABMD parameters to program units ---------------------
*=======================================================================

      IF(abmd) THEN
         CALL covfod(spring,abmd_tors)
      END IF

*=======================================================================
*----- Change the folding parameters to program units ------------------
*=======================================================================

      IF(lenerg) CALL covuen(uealfa,uemin,uemax)

*=======================================================================
*--------- Change switching function parameters to a suitable form -----
*=======================================================================

      hacut=(hanoff+hacut)*pi/180.0D0
      hanon=hanon*pi/180.0D0
      hanoff=hanoff*pi/180.0D0
      hana=hanoff
      hanb=hanon
      hanon=DCOS(hanon)
      hanoff=DCOS(hanoff)
      hacut=hana+hacut
      hacut=DCOS(hacut)
      hrcut=hrsoff+hrcut
      nhskip=nhskip+1

*=======================================================================
*---- Convert pressure variables ---------------------------------------
*=======================================================================

      CALL covpr(volume,t,taut,taup,pext,masspr,wpr,compressibility)

*================= END OF EXECUTABLE STATEMENTS ========================

1     FORMAT('=',78(' '),'='/80('=')/ /)
100   FORMAT(/'---> Testing Random Number: u =',e20.12,/
     &'---> It should be: u = 0.812053811384E-01')
      RETURN
      END
