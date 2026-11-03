      SUBROUTINE setup_box(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv
     &     ,xpa,ypa,zpa,xpb,ypb,zpb,index,cg,tg,list,hlist,rlist,nb,work
     &     ,iret,errmsg)

************************************************************************
*                                                                      *
*     Initialize the particle coordinates for the simulation.          *
*     For a simulation of a liquid, a first molecular configuration    *
*     is set up consistently with the chosen cubic lattice. The        *
*     molecular orientaion is not randomized. For a macromolecular     *
*     simulation, initial coordinates for the atoms are read in        *
*     free format. The routine checks for consistency with the         *
*     list of labels constructed in START.                             *
*     If the simulation involves both solvent molecules and            *
*     macromolecule, a macromolecule is inserted in an equilibrated    *
*     configuration of the liquid.                                     *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,nb
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),gh(*)
      CHARACTER*80 errmsg
      REAL*8  xpa(*),ypa(*),zpa(*),xpb(*),ypb(*),zpb(*),cg(3,3,*)
     &     ,tg(3,*),work(*)
      INTEGER mapnl(*),mapnl_slv(*),list(3,*),hlist(3,*),rlist(3,*)
     &     ,index(*)


*-------------------- LOCAL VARIABLES ----------------------------------

      EXTERNAL near0
      LOGICAL near0,mask(m1),mask2(m1)
      INTEGER i,j,n,nlist,tot,gnmol,pnmol,pnts,pntap,nmapnl,l1,offset
     &     ,nato_slt,tstep,tot1,tot2
      REAL*8 dummy,xpcm,ypcm,zpcm,xpcmp,ypcmp,zpcmp,fstep,cod(3,3)
      DATA dummy /0.0d0/

*==================== EXECUTABLE STATEMENTS ============================

************************************************************************
*                                                                      *
*    Set up box and topology if symmetry operation have be be          *
*    performed                                                         *
*                                                                      *
************************************************************************

      iret=0
      IF(nflag(1) .EQ. 0) THEN

*=======================================================================
*--------- Initialize coordinates for the Nose-Hover chain -------------
*=======================================================================

         IF(hoover) THEN
            DO i=1,nhoov
               gh(i)=0.0D0
            END DO
         END IF
      END IF

      IF(nflag(1) .NE. 0 .OR. analys) THEN
         IF(nflag(1) .EQ. 4) THEN
            DO i=1,nhoov
               gh(i)=0.0D0
            END DO
         END IF
         
*=======================================================================
*----- Read in template file -------------------------------------------
*=======================================================================

         IF(template) THEN
            nato_slt=ntap-nato_slv*nmol
            write(kprint,32000)
            CALL rdcmac(ktemplate,kprint,nres(1,1),beta,xpt0,ypt0
     &           ,zpt0,nato_slt,.FALSE.,1,1,iret,errmsg)
            IF(iret .EQ. 1) RETURN
         END IF

*=======================================================================
*----- Open direct access file for dumping -----------------------------
*=======================================================================

         IF(nconf.NE.0) THEN
            tstep=maxrun/nconf
         END IF

         IF(dmprnd) THEN
            n=ntap
            tot1=nbyte+9*2*nbyte
            tot2=n
            IF(n .LT. divide_records) THEN
               iret=1
               errmsg=
     &              'Cannot divide record by a number larger than'/ /
     &              ' the number of data. Abort.'
               RETURN
            END IF

            IF(analys) tstep=-1
            CALL open_file_dump(kcnfig,dmpfil,tot1,tot2,divide_records
     &           ,atom_record,nbyte,tstep,occupy_space,analys,iret
     &           ,errmsg)
            IF(iret .EQ. 1) RETURN
            IF(occupy_space .AND. (.NOT. analys)) THEN
               write(kprint,50000)
               DO i=1,3*ntap
                  work(i)=0.0D0
               END DO
               fstep=0.0D0
               DO i=1,3
                  DO j=1,3
                     cod(i,j)=0.0D0
                  END DO
               END DO
               DO i=1,tstep
                  CALL write_confc(cod,work(1),work(ntap+1),work(2*ntap
     &                 +1),ntap,fstep,i,1,divide_records,atom_record)
               END DO
            END IF
         END IF 
         RETURN
      END IF

*=======================================================================
*--------- Read group information from file sgroup.dat -----------------
*=======================================================================

      write(kprint,34000)

      IF(slv_add) THEN
         CALL add_solvent_tpg(mapnl,mapnl_slv,iret,errmsg)
         IF(iret .EQ. 1) RETURN
      END IF
      pntap=ntap
      IF(sgroup .AND. slt_exist) THEN
         CALL read_space_group(kprint,kgroup,cgroup,gnmol,cg,tg,iret
     &        ,errmsg)
         IF(iret .EQ. 1) RETURN

*=======================================================================
*--- Replicate topology of the structure -------------------------------
*=======================================================================

         CALL repl_tpg(mapnl,gnmol*icl*icm*icn,iret,errmsg)
         IF(iret .EQ. 1) RETURN
         write(kprint,30500)
         write(kprint,33000) ntap,lbond,lstretch,lconstr
     &        ,lbend,ltors,litor,int14p
         write(kprint,31500)
      END IF

*=======================================================================
*----- Write the topological file when required ------------------------
*=======================================================================

      IF(stoprun .AND. tpgwbn .AND. (slv_add .OR. sgroup)) THEN
         CALL set_ss_array(ss_point,m1+1,ss_index,nmol,nato_slv,ntap)
         CALL min_pack(ntap,mapnl,nmapnl)
         CALL dumptp(ktpgprm_write,mapnl,nmapnl)
      END IF
      IF(stoprun) RETURN

************************************************************************

*=======================================================================
*----- Read in the coordinates of the system if needed -----------------
*=======================================================================

      IF(pdb_read) THEN
         write(kprint,32300)
         CALL rdcmac(kconf,kprint,nres(1,1),beta,xp0,yp0,zp0
     &        ,pntap,.FALSE.,1,1,iret,errmsg)
         IF(iret .EQ. 1) RETURN
         
*=======================================================================
*----- Compute the positions of hydrogens if required ------------------
*=======================================================================

         CALL serchd(xp0,yp0,zp0,pntap,concta,m1,list,hlist,
     x        rlist,nlist,iret,errmsg)
         IF(iret.NE.0) RETURN
         CALL atthd(list,hlist,rlist,nlist,xp0,yp0,zp0,
     x        beta,pntap,xpa,ypa,zpa,iret,errmsg)
         IF(iret.NE.0) RETURN

      ELSE IF(slt_create) THEN

         write(kprint,32100)
         CALL rdcmac(kcoord_slt,kprint,nres(1,1),beta,xp0,yp0,zp0
     &        ,pntap,.FALSE.,1,1,iret,errmsg)
         IF(iret.NE.0) RETURN
         
*=======================================================================
*----- Compute the positions of hydrogens if required ------------------
*=======================================================================

         CALL serchd(xp0,yp0,zp0,pntap,concta,m1,list,hlist,
     x        rlist,nlist,iret,errmsg)
         IF(iret.NE.0) RETURN
         CALL atthd(list,hlist,rlist,nlist,xp0,yp0,zp0,
     x        beta,pntap,xpa,ypa,zpa,iret,errmsg)
         IF(iret.NE.0) RETURN
      END IF
         
*=======================================================================
*----- Read in template file -------------------------------------------
*=======================================================================

      IF(slt_exist) THEN
         IF(template) THEN
            write(kprint,32000)
            CALL rdcmac(ktemplate,kprint,nres(1,1),beta,xpt0,ypt0
     &           ,zpt0,pntap,.FALSE.,1,1,iret,errmsg)
            IF(iret.NE.0) RETURN
         END IF
      END IF
         
*=======================================================================
*--- Replicate coordinates ---------------------------------------------
*=======================================================================

      IF(slt_exist) THEN
         
         IF(sgroup) THEN
            CALL change_frame(co,oc,-1,pntap,xp0,yp0,zp0,xpa,ypa,zpa)
            CALL repl_coord(xpa,ypa,zpa,nb,pntap,cg,tg,gnmol,icl
     &           ,icm,icn,pnmol,pnts,iret,errmsg)
            IF(iret .EQ. 1) RETURN
            CALL change_frame(co,oc,1,pnts,xpa,ypa,zpa,xp0,yp0
     &           ,zp0)
            
*=======================================================================
*--- Replicate coordinates of the template structure -------------------
*=======================================================================

            IF(template) THEN
               CALL change_frame(co,oc,-1,pntap,xpt0,ypt0,zpt0,xpa,ypa
     &              ,zpa)
               CALL repl_coord(xpa,ypa,zpa,nb,pntap,cg,tg,gnmol
     &              ,icl,icm,icn,pnmol,pnts,iret,errmsg)
               IF(iret .EQ. 1) RETURN
               CALL change_frame(co,oc,1,pnts,xpa,ypa,zpa,xpt0
     &              ,ypt0,zpt0)
            END IF
         END IF
         CALL clpcm(mass,wmtp,xp0,yp0,zp0,ntap,1,xpcmp,ypcmp,zpcmp)
         
         write(kprint,10000) xpcmp,ypcmp,zpcmp
****   to be  remove         

*         CALL plotc(fstep,beta,co,xp0,yp0,zp0,ntap,nres,m1,prsymb)
*         IF(iret .EQ. 0) STOP

!=======================================================================
!------ Reset the position of the centre of mass to the centre of ------
!------ the simulation box and check the box size against the ----------
!------ dimension of the macromolecule ---------------------------------
!=======================================================================

         IF(rescm) THEN
            DO i=1,ntap
               xp0(i)=xp0(i)-xpcmp
               yp0(i)=yp0(i)-ypcmp
               zp0(i)=zp0(i)-zpcmp
            END DO
         END IF
      END IF

      

      IF(slv_create) THEN

*=======================================================================
*--------- Add solvent topology to the topology file -------------------
*=======================================================================

         IF(slv_generate) THEN

*=======================================================================
*----- Read in the solvent coordinates in pdb form ---------------------
*=======================================================================

            write(kprint,32200)
            CALL rdcmac(kcoord_slv,kprint,nres_slv(1,1),beta_slv
     &           ,xpa,ypa,zpa,nato_slv,.FALSE.,1,1,iret,errmsg)
            IF(iret.NE.0) RETURN
         
*=======================================================================
*----- Compute the positions of hydrogens if required ------------------
*=======================================================================

            CALL serchd(xpa,ypa,zpa,nato_slv,concta_slv,slvatm,list
     &           ,hlist,rlist,nlist,iret,errmsg)
            IF(iret.NE.0) RETURN

            CALL atthd(list,hlist,rlist,nlist,xpa,ypa,zpa,
     &           beta_slv,nato_slv,xpb,ypb,zpb,iret,errmsg)
            IF(iret.NE.0) RETURN
            CALL generate_slv(co,oc,xpa,ypa,zpa,xpb,ypb,zpb,mass_slv
     &           ,nato_slv,nmol,icl_slv,icm_slv,icn_slv,nform
     &           ,rmol,boxl,slv_randomize,iret,errmsg)
            IF(iret.NE.0) RETURN

            IF(linser) THEN
               CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xp0,yp0,zp0)
               CALL change_frame(co,oc,-1,nmol*nato_slv,xpa,ypa,zpa,xpa
     &              ,ypa,zpa)
               CALL insert_solute(kprint,co,xp0,yp0,zp0,xpa,ypa
     &              ,zpa,xpb,ypb,zpb,index,nbtype,radius,pnbd1
     &              ,nbtype_slv,ntap,nmol,nato_slv)
               CALL change_frame(co,oc,1,ntap,xp0,yp0,zp0,xp0,yp0,zp0)
               CALL change_frame(co,oc,1,nmol*nato_slv,xpa,ypa,zpa,xpa
     &              ,ypa,zpa)
            END IF
         END IF
         IF(slv_read) THEN
            offset=1
            write(kprint,32200)
            DO i=1,nmol
              CALL rdcmac(kcoord_slv,kprint,nres_slv(1,1)
     &             ,beta_slv,xpa(offset),ypa(offset),zpa(offset)
     &             ,nato_slv,.TRUE.,i,nmol,iret,errmsg)
               IF(iret.NE.0) RETURN
               CALL serchd(xpa(offset),ypa(offset),zpa(offset)
     &              ,nato_slv,concta_slv,slvatm,list,hlist,rlist,nlist
     &              ,iret,errmsg)
               IF(iret.NE.0) RETURN
               CALL atthd(list,hlist,rlist,nlist,xpa(offset)
     &              ,ypa(offset),zpa(offset),beta_slv,nato_slv,xpb,ypb
     &              ,zpb,iret,errmsg)
               IF(iret.NE.0) RETURN
               offset=offset+nato_slv
            END DO
            IF(linser) THEN
               
               CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xp0,yp0,zp0)
               CALL change_frame(co,oc,-1,nmol*nato_slv,xpa,ypa,zpa,xpa
     &              ,ypa,zpa)
               CALL insert_solute(kprint,co,xp0,yp0,zp0,xpa,ypa
     &              ,zpa,xpb,ypb,zpb,index,nbtype,radius,pnbd1
     &              ,nbtype_slv,ntap,nmol,nato_slv)
               CALL change_frame(co,oc,1,ntap,xp0,yp0,zp0,xp0,yp0,zp0)
               CALL change_frame(co,oc,1,nmol*nato_slv,xpa,ypa,zpa,xpa
     &              ,ypa,zpa)
            END IF
         END IF
         CALL add_solvent_coord(xp0,yp0,zp0,xpa,ypa,zpa,ntap,nmol
     &        ,nato_slv,m1,iret,errmsg)
         IF(iret.NE.0) RETURN

         CALL add_solvent_tpg(mapnl,mapnl_slv,iret,errmsg)
         IF(iret.NE.0) RETURN

      END IF


      IF(slt_create .OR. slv_create .OR. slv_add) THEN

         CALL set_ss_array(ss_point,m1+1,ss_index,nmol,nato_slv,ntap)

         IF(tpgwbn) THEN
            CALL min_pack(ntap,mapnl,nmapnl)
            CALL dumptp(ktpgprm_write,mapnl,nmapnl)
         END IF

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
         IF(iret.NE.0) RETURN

*=======================================================================
*----- Compute the total charge and scale it to zero when required -----
*=======================================================================

         CALL scale_charges(kprint,nprot_charges,prot_charges
     &        ,chrge,protl,nprot,scharge)

*=======================================================================
*----- Change charges to program units ---------------------------------
*=======================================================================

         CALL covchg(chrge,ntap)
               
*=======================================================================
*--------- Change the covalent interaction parameters ------------------
*---------------- program units ----------------------------------------
*=======================================================================
            
         CALL covbd(potbo,potbe,potto,potit,ptorj,lstretch,lbend
     &        ,ltors,litor,itor_ptype,m9,m2,m3,m4)

*=======================================================================
*--------- Change the non bonded interaction parameters to program -----
*--------- units and compute the macromolecule total mass --------------
*=======================================================================

         CALL covnb(pnbd1,pnbd2,pnbd3,pnbd4,sjorg,ejorg,lpnbd,iz
     &        ,phyd1,phyd2,lphyd,mass,wmtp,ntap,ecc6,ecc12,ecc146
     &        ,ecc1412,c6jorg,c12jorg,type_table,lj_fudge,m6)
               
*=======================================================================
*--------- Compute reduced mass for atomic groups ----------------------
*=======================================================================

         CALL redmss(ngrp,grppt,mass,pmass)

*=======================================================================
*--------- Compute the intergroup interaction map ----------------------
*=======================================================================
         
         CALL igmap(ngrp,grppt,ingrpp,ingrp,m12,mapnl,errmsg,iret)
         IF(iret.EQ.1) RETURN
            
*=======================================================================
*--------- Compute pointers to start an end of residue -----------------
*=======================================================================
         
         CALL pntres(nres(1,1),mres,atomg,grppt,nbun,resg,ngrp,ntap)
      END IF

      write(kprint,40000) nmol,nato_slv,nprot-nmol,ntap
     &     ,lbond,lstretch,lconstr,lbend,ltors,litor,int14p

*=======================================================================
*----- Open direct access file for dumping -----------------------------
*=======================================================================

      IF(nconf.NE.0) THEN
         tstep=maxrun/nconf
      END IF
      IF(dmprnd) THEN
         IF(analys) tstep=-1
         n=ntap
         tot1=nbyte+9*2*nbyte
         tot2=n
         IF(n .LT. divide_records) THEN
            iret=1
            errmsg=
     &           'Cannot divide record by a number larger than'/ /
     &           ' the number of data. Abort.'
            RETURN
         END IF

         CALL open_file_dump(kcnfig,dmpfil,tot1,tot2,divide_records
     &        ,atom_record,nbyte,tstep,occupy_space,analys,iret,errmsg)
         IF(iret .EQ. 1) RETURN
         IF(occupy_space .AND. (.NOT. analys)) THEN
            write(kprint,50000)
            DO i=1,3*ntap
               work(i)=0.0D0
            END DO
            fstep=0.0D0
            DO i=1,3
               DO j=1,3
                  cod(i,j)=0.0D0
               END DO
            END DO
            DO i=1,tstep
               CALL write_confc(cod,work(1),work(ntap+1),work(2*ntap+1)
     &              ,ntap,fstep,i,1,divide_records,atom_record)
            END DO
         END IF
      END IF 
      
*================= END OF EXECUTABLE STATEMENTS ========================

10000 FORMAT(
     &     /'                    Solute Center of Mass coordinates at '
     &     /'                    ',3f12.3/)
20000 FORMAT(//'Protein Center of Mass Reset to zero --->'//)
30500 FORMAT(//
     &'                    ------------------------------------------'
     &     //
     &'                                In The New System             ')
31000 FORMAT(//'                Found --->',i6,' Atoms     ',
     &     i6,' Bonds     ',i6,' Angles    '/
     &     '                          ',i6,' P-Torsions',
     &     i6,' I-Torsions',i6,' 1-4 Interactions')
31500 FORMAT(//
     &'                    ------------------------------------------')
32000 FORMAT(//
     &'                    =========================================='/
     &'                    =                                        ='/
     &'                    =            Read Template File          ='/
     &'                    =            to Compute X-rms            ='/
     &'                    =                                        ='/
     &'                    =========================================='
     &     //)
321000 FORMAT(//
     &'                    =========================================='/
     &'                    =        Reading Solute Coordinates      ='/
     &'                    =========================================='
     &     //)
322000 FORMAT(//
     &'                    =========================================='/
     &'                    =       Reading Solvent Coordinates      ='/
     &'                    =========================================='
     &     //)
323000 FORMAT(//
     &'                    =========================================='/
     &'                    =       Reading System Coordinates       ='/
     &'                    =========================================='
     &     //)
33000 FORMAT(//
     &  '  *****************************************************',
     &  '*******************'/
     &  '  *                          TOPOLOGY  List            ',
     &  '                  *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *          ',i6,' Atoms        ',i6,' Bonds     ',i6,
     &  ' FLexible Bonds  *'/
     &  '  *          ',i6,' Rigid Bonds  ',i6,' Angles    ',i6,
     &  ' P-Torsions      *'/
     &  '  *          ',i6,' I-Torsions   ',i6,' 1-4 Inter.',6x,
     &  '                 *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *****************************************************',
     &  '*******************'/)
34000 FORMAT(/
     &  '   =======================================================',
     &  '================='/
     &  '   =                                                      ',
     &  '                ='/
     &  '   =                 Setting up the Simulation Box        ',
     &  '                ='/
     &  '   =                                                      ',
     &  '                ='/
     &  '   =======================================================',
     &  '================='/)
40000 FORMAT(//
     &  '  *****************************************************',
     &  '*******************'/
     &  '  *                        New TOPOLOGY  List          ',
     &  '                  *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *          ',i6,' Solvent Mol  ',i6,' Atoms Each',i6,
     &  ' Solute Mol.     *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *                            For the Sysyem          ',
     &  '                  *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *          ',i6,' Atoms        ',i6,' Bonds     ',i6,
     &  ' FLexible Bonds  *'/
     &  '  *          ',i6,' Rigid Bonds  ',i6,' Angles    ',i6,
     &  ' P-Torsions      *'/
     &  '  *          ',i6,' I-Torsions   ',i6,' 1-4 Inter.',6x,
     &  '                 *'/
     &  '  *                                                    ',
     &  '                  *'/
     &  '  *****************************************************',
     &  '*******************'/)
50000 FORMAT(/'        <------ Zeroing Trajectory File ------->'/
     &        '              This might take a while ...       '/)

      RETURN
      END

      SUBROUTINE int_str (icurr,string,n)
c     get integer ICURR of length n and return string STRING of
c     same length; icurr is damaged on return
      
      INTEGER icurr,ibase,i,idum(8),ibase_max
      CHARACTER*8 string
      CHARACTER*11 ints
      ints="0123456789"
      string="        "
      do n=0,7
         ibase_max = icurr/10**n 
         if(ibase_max.eq.0) THEN
            go to 1
         END IF
      end do
 1    do i=1,n
         ibase = 10**(n-i)
         idum(i) = icurr/ibase
         icurr = icurr - ibase*idum(i)
      end do
      
      do i=1,n
         ibase=idum(i)+1
         string(i:i) = ints(ibase:ibase)
      end do

      RETURN
      END 
