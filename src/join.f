      SUBROUTINE join(nbun,ntap,ngrp,lbond,lbend,ltors,litor,lphyd,lpnbd
     &     ,nbone,int14p,int13p,llacc,lldon,mback,grppt,int14,int13,nres
     &     ,mend,lacc,ldon,nhtype,lbnd,lbndg,ltor,litr,atres,concta
     &     ,nbtype,mass,beta,betb,alnbd,potbo,potbe,potto,potit,chrge
     &     ,bsitp,asitp,nbsitp,nasitp,nrigg,prsymb,mm1,mm2,mm3,mm4,mm6
     &     ,mm9,mm10,mm11,mm13,mm14,debug_bt,debug_pt,debug_it,debug_rs
     &     ,debug_ct,debug_st,adihed,mesg)

************************************************************************
*   Time-stamp: <97/02/08 13:42:37 marchi>                             *
*                                                                      *
*   mm1:  Dimension  of the list of atoms                              *
*   mm2:  Dimension of the list of bendings                            *
*   mm3:  Dimension of the list of torsions                            *
*   mm4:  Dimension of the list of improper torsions                   *
*   mm6:  Dimension of the types of hydrogen bond (obsolescent)        *
*   mm9:  Dimension of the list of bonds                               *
*   mm10: 2nd dimension of the connection list (i.e. how many          *
*         connection per atom are allowed)                             *
*   mm11: Dimension of the list of groups                              *
*   mm13: 1st dimension of bsitep and asitep                           *
*   mm14: 2nd dimension of bsitep and asitep                           *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Feb  7 1997 -                                     *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER mm1,mm2,mm3,mm4,mm6,mm9,mm10,mm11,mm13,mm14
      INTEGER nbun,ntap,ngrp,lbond,lbend,ltors,litor,llacc,lldon,lphyd
     &     ,lpnbd,int14p,int13p,nbone,mback(*),int14(2,*),int13(2,*)
     &     ,grppt(2,*),nhtype(mm6,*),nres(mm1,*),mend(*),lbnd(2,*)
     &     ,lbndg(3,*),ltor(4,*),litr(4,*),concta(mm1,*),nbtype(*)
     &     ,lacc(2,*),ldon(2,*),bsitp(mm13,*),asitp(mm13,*),nbsitp(*)
     &     ,nasitp(*),atres(2,*),nrigg,nrigat
      REAL*8 mass(*),potbo(mm9,*),potbe(mm2,*),potto(mm3,*),potit(mm4,*)
     &     ,chrge(*)
      CHARACTER*7 beta(*),betb(*),alnbd(*)
      CHARACTER*8 prsymb(*)
      LOGICAL debug_bt,debug_pt,debug_it,debug_rs,debug_ct,debug_st
     &     ,adihed


*-------------------- PARAMETER STATEMENT ------------------------------

      INCLUDE 'parameters.h'

      INTEGER iret,bugs(1001),nbugs,b1,b2,b3,b4
      logical print_bond(1001) 
      character*28 types1,types2
      CHARACTER*80 errmsg
      INTEGER lben(3,m1),lbens      
      CHARACTER*7 clinks(nores,2)
      CHARACTER*9 mesg

*-----------------------------------------------------------------------

      INTEGER i,ia,ib,j,j1,k,m,l,n,npa,npan,nato,nta,ntd,nme,
     x        nbonds,nlinks,ntph,connct(noato,m10)
      INTEGER coordi,ncount,nchya,nchyd,lconstr,notfound
      REAL*8  sum

      INTEGER protl(2*m1),nprot
      LOGICAL mask(m1),mask2(m1)
      COMMON /rag1/ nprot,protl,mask,mask2
*==================== EXECUTABLE STATEMENTS ============================

      notfound=0
      iret=0

*=======================================================================
*----- Welcome message to the routine ----------------------------------
*=======================================================================

      WRITE(kprint,20000) mesg

*=======================================================================
*----- Set up an integer table for all the hydrogen bond interactions
*=======================================================================

      CALL phapp(alhyd,MM6,lphyd,alnbd,lpnbd,nhtype,iret,errmsg)

*----------------- Iret is 1 STOP the run ------------------------------

      IF(iret.EQ.1) THEN
          WRITE(kprint,'(a)') errmsg
          STOP
      END IF
      IF(nbun .EQ. 0) RETURN

*=======================================================================
*------------ Begin to loop over the solute residues -------------------
*=======================================================================

      ncount=0
      nchya=0
      nchyd=0
      nbone=0
      ngrp=0
      DO i=1,nbun
         npa=mend(i)
         IF(i.NE.nbun) THEN
            npan=mend(i+1)
         ELSE
            npan=0
         END IF
         nato=natop(npa)
         nbonds=nmbo(npa)
         atres(1,i)=ncount+1
         atres(2,i)=ncount+nato

*----- Hydrogen bonds

         nta=nacc(npa)
         ntd=ndon(npa)

*=======================================================================
*----------------- Append acceptor names to the list -------------------
*=======================================================================

         CALL hapacp(alpha(1,npa),NOATO,nato,jacc(1,1,npa),HYDP,
     x        lacc,MM1,nta,nchya,ncount,iret,errmsg)

*----------------- If iret .EQ. 1 STOP !! ------------------------------
*                                 ====

         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*-----------------------------------------------------------------------


*=======================================================================
*---------------- Append donor names to the list -----------------------
*=======================================================================
         
         CALL hapdon(alpha(1,npa),NOATO,nato,jdon(1,1,npa),HYDP,
     x                ldon,MM1,ntd,nchyd,ncount,iret,errmsg)

*----------------- If iret .EQ. 1 STOP !! ------------------------------
*                                ====

         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*-----------------------------------------------------------------------


*=======================================================================
*------ Append lists of the atom labels of the residue to the
*------ protein lists
*=======================================================================


         CALL appatm(alpha(1,npa),natype(1,npa),qchge(1,npa),qmss(1,npa)
     &        ,npa,i,nato,beta,betb,chrge,mass,nres,MM1,ncount)

*=======================================================================
*------ Create a list of pointers to charge groups ---------------------
*=======================================================================

         CALL appgrp(ngrp,grppt,MM11,jngrp(npa),jgrppt(1,1,npa),ncount,
     x         iret,errmsg)

*----------------- If iret .EQ. 1 STOP !! ------------------------------
*                                 ====

         IF(iret.EQ.1) THEN
            CALL xerror(errmsg,80,1,2)
         END IF
*-----------------------------------------------------------------------

*=======================================================================
*------ Append list of backbone atoms of the residue to the ------------
*------ protein lists --------------------------------------------------
*=======================================================================

         nme=nback(npa)
         CALL appbkn(alpha(1,npa),nato,alphb(1,npa),nme,mback,ncount,
     x        nbone,iret,errmsg)
         
*----------------- If iret .EQ. 1 STOP ** ------------------------------
*                                 ====

         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
*-----------------------------------------------------------------------


*=======================================================================
*------ Prepare the connection table for the residue and append
*------ it to the one of the protein. No connection between units
*------ are included at this stage.
*=======================================================================

         CALL sctab(alpha(1,npa),nato,jbnd(1,1,npa),NOATO,nbonds,
     x        connct,NOATO,MM10)
         CALL appct(connct,concta,NOATO,MM1,MM10,nato,ncount)

*-----------------------------------------------------------------------

*=======================================================================
*----- Create a table which contains the links between each two
*----- residues.
*=======================================================================

         IF(i.NE.nbun) THEN
            IF(i.EQ.1) THEN
               clinks(i,1)=cnat(2,npa)
               clinks(i,2)=cnat(1,npan)
            ELSE
               clinks(i,1)=cnat(2,npa)
               clinks(i,2)=cnat(1,npan)
            END IF
         END IF

*------ Increase counters


         nbone=nbone+nme
         ncount=ncount+nato
         nchya=nchya+nta
         nchyd=nchyd+ntd
         IF(ncount .GT. sitslu) THEN
            errmsg=' Number of atoms exceed physical '//
     x           ' dimensions: increase _SIT_SOLU_ in config.H'
            iret=1
            CALL xerror(errmsg,80,1,2)
         END IF

      END DO

      IF(debug_rs) THEN
        WRITE(kprint,*)
        WRITE(kprint,'(a)') '============ RESIDUES SEQUENCE ======'//
     x       '======'
        WRITE(kprint,'(a8,2x,i6)') (prsymb(mend(i)),i,i=1,nbun)
      END IF


*----- Number of atoms, acceptors, donors and links in the protein

      llacc=nchya
      lldon=nchyd
      ntap=ncount
      nlinks=nbun-1

*=======================================================================
*-------------- Link the protein together ------------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Linking the molecule together      ---->'
      CALL linka(beta,ncount,clinks,nores,nlinks,concta,MM1,nres(1,1))

*=======================================================================
*----- Add extra links -------------------------------------------------
*=======================================================================

      IF(mesg .EQ. 'SOLUTE:  ') THEN
         WRITE(kprint,'(5x,22a)')
     &        mesg,'Adding extra links                 ---->'
         CALL addlnk(beta,atres,concta,MM1,MM10,rbond,xbond,xnbond,
     x        iret,errmsg)
         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*---- Find all possible bonds
*=======================================================================

      WRITE(kprint,'(5x,22a)')
     &     mesg,'Searching for Bonds                ---->'
      CALL sbond(concta,MM1,ncount,lbnd,lbond,MM9,iret)

*====== DEBUG WRITE statements =========================================

      IF(debug_ct) THEN
      WRITE(kprint,*)
      WRITE(kprint,'(a)')
     &     '============ ATOMS of the macromolecule ========'
      sum=0.0d0
      DO 120 i=1,ntap
        WRITE(kprint
     &       ,'(1h ,i4,2x,a4,2x,a4,2x,f7.1,2x,i4,f7.3,5x,i4,4x,i4)' )i
     &       ,beta(i),betb(i),mass(i),nbtype(i),chrge(i),nres(i,1)
     &       ,nres(i,2)
          sum=sum+chrge(i)
120   CONTINUE
      WRITE(kprint,'(20x,18h   Total charge = ,f10.4)')
     &     sum
      WRITE(kprint,'(5x,a)')'=== CONNECTION TABLE ==='

          DO 131 i=1,ntap
              coordi = concta(i,1)
              WRITE(kprint,'(1h ,a4,i6,a,i1,a,5(i6
     &,2x))')beta(i),i,' (',coordi,') : ',(concta(i,j),j=2,coordi+1)
131       CONTINUE

      WRITE(kprint,*)
      WRITE(kprint,'(5x,a)'
     &     )'======= BONDS TABLE ======= '

      DO 140 m=1,lbond
        WRITE(kprint,'(1h ,a4,2x,a4,5x,i6,2x,i6,10x
     &,i6)')beta(lbnd(1,m)),beta(lbnd(2,m)),lbnd(1,m),lbnd(2,m),m
140   CONTINUE
      END IF
      IF (iret .EQ. 1) STOP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



*=======================================================================
*---- Find all possible bendings ---------------------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Searching for Bendings             ---->'
      CALL sbend(concta,MM1,ncount,lbndg,lbend,MM2)

*=======================================================================
*---- Find angle bend to be excluded  ----------------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Searching for Bendings to Omit     ---->'
      CALL pckben(nbun,mend,natop,beta,jbend,nbend,NOATO,lben,
     x            lbens,iret,errmsg)

      IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*=======================================================================
*---- Assign a type to each atom of the protein ------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Assigning a Type to Each Atom      ---->'
      CALL asstpe(betb,alnbd,massu,ntap,lpnbd,nbtype,mass,bugs,nbugs)

      if(nbugs.ne.0)  THEN 
         call xerror(errmsg,80,1,221)
         errmsg = ' In The TPG file the following atoms have a type'
     &        //  ' which'
         call xerror(errmsg,80,1,222)
         errmsg = ' was not found among those given in the PRM file'
         call xerror(errmsg,80,1,222)
         errmsg = '  atom |  label   | type     |  residue'
         call xerror(errmsg,80,1,222)
         do ia = 1,nbugs
            print_bond(ia)=.true.
         end do
         DO ia = 1,nbugs
            k=bugs(ia) 
c---        see if type has been printed already
            if(print_bond(ia)) THEN 
               do ib =ia,nbugs
                  m=bugs(ib) 
                  if(betb(k).eq.betb(m)) print_bond(ib) = .false.
               end do
               b1=mend(nres(k,1)) 
               write(kprint,1007) k,beta(k),betb(k)
     &              ,nres(k,1),prsymb(b1)
 1007          format(i8,2x,'|',2(2x,a4,4x,'|'),1x,i4,a4)
            END IF
         END DO 
         IF(nbugs.eq.1001) THEN
            errmsg = ' ..... and even more!!' 
            call xerror(errmsg,80,1,222)
         END IF
         call xerror(errmsg,80,1,223)
         STOP
      END IF
      if(debug_ct) THEN 
        WRITE(kprint,*)
        WRITE(kprint,1035)
1035    FORMAT('============ ATOMS of the macromolecule ====='
     &       /' =====(after mass assigmente from prm file  =====')
        DO i=1,ntap
          WRITE(kprint
     &         ,'(1h ,i4,2x,a4,2x,a4,2x,f7.1,2x,i4,f7.3,5x,i4,4x,i4)' )i
     &         ,beta(i),betb(i),mass(i),nbtype(i),chrge(i),nres(i,1)
     &         ,nres(i,2)
          sum=sum+chrge(i)
        END DO
      END IF

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


*=======================================================================
*---- Find all possible 1-4 interactions -------------------------------
*---- Copy the 1-4 interactions to an array ----------------------------
*=======================================================================

      CALL chnpr(protl,2*m1,npm,nprot,concta,mm1,mask,mask2,ntap
     &     ,iret,errmsg)
      IF(iret .EQ. 1) call xerror(errmsg,80,1,222)

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Find All Possible 1-4 Interactions ---->'
      CALL clt14(concta,MM1,MM10,ntap,lbndg,lbend,int14,int14p,MM3,ltor,
     x           ltors,nprot,protl,adihed)

*=======================================================================
*----- Produce a list of 1-3 interactions ------------------------------
*=======================================================================

      int13p=lbend
      DO i=1,lbend
          int13(1,i)=lbndg(1,i)
          int13(2,i)=lbndg(3,i)
      END DO

*=======================================================================
*---- Remove omitted angle bend from the original list -----------------
*=======================================================================

      CALL rmbend(lbndg,lbend,lben,lbens)



*=======================================================================
*---- Copy two 1-dimensional parameter arrays to -----------------------
*---- a 2-dimensional array --------------------------------------------
*=======================================================================

      DO j1=1,lpbon
          pbon(j1,1)=pbon1(j1)
          pbon(j1,2)=pbon2(j1)
      END DO

*=======================================================================
*----  Match the bond potential parameters -----------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Searching for Bond Parameters  ---->'
      CALL pbond(lbnd,lbond,MM9,betb,albon,pbon,lpbon,BONDP,
     x     potbo,2,bugs,nbugs)

      IF(debug_st) THEN
         WRITE(kprint,*)
         WRITE(kprint,'(5x,a)'
     &        )'======= BONDS TABLE ======= '
    
         DO m=1,lbond
            WRITE(kprint
     &           ,'(1h ,a4,2x,a4,5x,i6,2x,i6,2x,f7.2,f10.3,10x,i4)'
     &           )betb(lbnd(1,m)),betb(lbnd(2,m)),lbnd(1,m),lbnd(2,m)
     &           ,potbo(m,1),potbo(m,2),m
         END DO
      END IF

*----------------- If iret .EQ. 1 STOP ** ------------------------------

      if(nbugs.ne.0)  THEN 
         call xerror(errmsg,80,1,221)
         errmsg = 'Bond parameters not found for the following pair:'
         call xerror(errmsg,80,1,222)
         errmsg = '            pair           |  types   | labels'//
     &        '   |        res.'
         call xerror(errmsg,80,1,222)
         do ia = 1,nbugs
            print_bond(ia)=.true.
         end do
c====    Start printing unassigned bonds 
         DO ia = 1,nbugs
            iret=bugs(ia) 
            k=lbnd(1,iret)
            m=lbnd(2,iret)
            types1 = betb(k)//betb(m)
c---        see if bond has been printed 
            if(print_bond(ia)) THEN 
               do ib = ia , nbugs
                  iret=bugs(ib) 
                  b1=lbnd(1,iret)
                  b2=lbnd(2,iret)
                  types2 = betb(b1)//betb(b2)
                  if(types1.eq.types2) print_bond(ib) = .false.
               end do
               b1=mend(nres(k,1)) 
               b2=mend(nres(m,1)) 
               write(kprint,1008) k,m, betb(k),betb(m),
     &              beta(k),beta(m),nres(k,1),prsymb(b1),nres(m,1)
     &              ,prsymb(b2)   
1008           format(" *BND",2i9,6x,' |',2(1x,a4),'|',2(1x,a4),'|',2(i4
     &              ,a3))      
            END IF
         END DO 
         IF(nbugs.eq.1001) THEN 
            errmsg = ' ..... and even more!!' 
            call xerror(errmsg,80,1,222)
         END IF
         call xerror(errmsg,80,1,223)
         notfound=1
      END IF


*----------------- If iret .EQ. 1 STOP ** ------------------------------

   

*------- Copy 4 1-dimensional parameter arrays to ----------------------
*---------------- a 4-dimensional array --------------------------------

      DO 60 j1=1,lpbnd
          pbnd(j1,1)=pbnd1(j1)
          pbnd(j1,2)=pbnd2(j1)
          pbnd(j1,3)=pbnd3(j1)
          pbnd(j1,4)=pbnd4(j1)
          pbnd(j1,5)=pbnd5(j1)
          pbnd(j1,6)=pbnd6(j1)
          pbnd(j1,7)=pbnd7(j1)
60    CONTINUE

*-------------- Match the bending potential parameters -----------------
*----------- Iret is equal to 1 if the match is incomplete -------------

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Searching for Bendings Parameters  ---->'
      CALL pbend(lbndg,lbend,MM2,betb,albnd,pbnd,lpbnd,BENDP,
     x           potbe,7,bugs,nbugs)

      IF(debug_bt) THEN
          WRITE(kprint,*)
          WRITE(kprint,'(a)'
     &         )'======= BENDINGS TABLE ======= '

          DO 70 m=1,lbend
            WRITE(kprint,'(1h ,a4,2x,a4,2x,a4,5x,i6,2x,i6,2x,i6,
     x                  2x,f7.2,f10.3,10x,i4)')
     x              betb(lbndg(1,m)),betb(lbndg(2,m)),betb(lbndg(3,m)),
     x              lbndg(1,m),lbndg(2,m),lbndg(3,m),potbe(m,1),
     x              potbe(m,2),m
70        CONTINUE
      END IF

      if(nbugs.ne.0)  THEN 
         call xerror(errmsg,80,1,221)
         errmsg = 'Bend parameters not found for the'// 
     &        ' following triplets - '
         call xerror(errmsg,80,1,222)
         errmsg = '           triplet         |    types      |'//
     &        '   labels      |       residues  '
         call xerror(errmsg,80,1,222)
         do ia = 1,nbugs
            print_bond(ia)=.true.
         end do
c====    Start printing unassigned bonds 
         DO ia = 1,nbugs
            iret=bugs(ia) 
            k=lbndg(1,iret)
            m=lbndg(2,iret)
            l=lbndg(3,iret)
            types1 = betb(k)//betb(m)//betb(l)
c---        see if bond has been printed 
            if(print_bond(ia)) THEN 
               do ib = ia , nbugs
                  iret=bugs(ib) 
                  b1=lbndg(1,iret)
                  b2=lbndg(2,iret)
                  b3=lbndg(3,iret)
                  types2 = betb(b1)//betb(b2)//betb(b3) 
                  if(types1.eq.types2) print_bond(ib) = .false.
               end do
               b1=mend(nres(k,1)) 
               b2=mend(nres(l,1)) 
               b3=mend(nres(m,1)) 
               write(kprint,1009) k,m,l, betb(k),betb(m)
     &              ,betb(l),beta(k),beta(m),beta(l),nres(k,1)
     &              ,prsymb(b1),nres(m,1),prsymb(b2),nres(l,1)
     &              ,prsymb(b3)
1009           format(" *ANG",3i7,3x,' |',3(1x,a4),'|',3(1x,a4),'|',3(i4
     &              ,a3))      
            END IF
         END DO 
         IF(nbugs.eq.1001) THEN 
            errmsg = ' ..... and even more!!' 
            call xerror(errmsg,80,1,222)
         END IF
         call xerror(errmsg,80,1,223)
         notfound=1
      END IF

*=======================================================================
*---- Pick torsions that match input -----------------------------------
*=======================================================================

      IF(.NOT. adihed) THEN
          WRITE(kprint,'(5x,2a)')
     &
     &        mesg,'Pick P-Torsions that Match Input   ---->'
          CALL pcktor(nbun,mend,natop,beta,jtor,ntor,NOATO,ltor,ltors,
     x            iret,errmsg)
          IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*---- Pick the rigid fragments that match input ------------------------
*=======================================================================

      WRITE(kprint,'(5x,22a)')
     &     mesg,'Pick rigid fragments that Match Input  ---->'

      CALL pckfrg(nbun,mend,natop,beta,stpr,stsc,nstpr,nstsc,nstrg,
     x     norsre,noatre,bsitp,asitp,nbsitp,nasitp,nrigg,MM13,
     x     MM14,iret,errmsg)
      nrigat=0
      DO i=1,nrigg
          nrigat=nrigat+nbsitp(i)
      END DO

      IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*---- Add extra Torsions -----------------------------------------------
*=======================================================================

      IF(mesg .EQ. 'SOLUTE:  ') THEN
         WRITE(kprint,'(5x,22a)')
     &        mesg,'Add extra P-Torsions               ---->'
         CALL addtor(beta,atres,rtor,xtor,xntor,ltor,ltors,MM3,
     x        .FALSE.,concta,MM1,iret,errmsg)
         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*------- Copy the two 1-dimensional parameter arrays to ----------------
*---------------- a 2-dimensional array --------------------------------
*=======================================================================

      DO 90 j1=1,lptor
          ptor(j1,1)=ptor1(j1)
          ptor(j1,2)=DFLOAT(ntor2(j1))
90    CONTINUE

*=======================================================================
*-------------- Match the torsion potential parameters -----------------
*----------- Iret is equal to 1 if the match is incomplete -------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Match P-Torsions with Parameters   ---->'

      CALL pptors(ltor,ltors,MM3,betb,altor,ptor,lptor,TORSP,
     x     potto,2,bugs,nbugs)


*----------------- If iret .EQ. 2 a WARNING is provided ----------------

      IF(debug_pt) THEN
         WRITE(kprint,*)
         WRITE(kprint,'(a)'
     &        )'======= PROPER TORSION TABLE ======= '
         
         DO 160 m=1,ltors
            ntph=DINT(potto(m,2)+0.5d0)
            WRITE(kprint
     &           ,'(1h ,4(a4,2x),2x,4(i4,1x),f7.2,i5,5x,i4)'
     &           )betb(ltor(1,m)),betb(ltor(2,m)),betb(ltor(3,m))
     &           ,betb(ltor(4,m)),ltor(1,m),ltor(2,m),ltor(3,m),ltor(4
     &           ,m),potto(m,1),ntph,m
 160     CONTINUE
      END IF
      
c==== If unassigned torsion are found print them out and stops 

      if(nbugs.ne.0)  THEN 
         call xerror(errmsg,80,1,221)
         errmsg = 'Torsion parameters not found for the'// 
     &        ' following quadruplets  - '
         call xerror(errmsg,80,1,222)
         errmsg = '               quadruplet        |'
     &        //'       types        |     labels         |  res     '
         call xerror(errmsg,80,1,222)
         do ia = 1,nbugs
            print_bond(ia)=.true.
         end do
c====    Start printing unassigned bonds 
         DO ia = 1,nbugs
            iret=bugs(ia) 
            k=ltor(1,iret)
            m=ltor(2,iret)
            l=ltor(3,iret)
            n=ltor(4,iret)
            types1 = betb(k)//betb(m)//betb(l)//betb(n)
c---        print torsion ia and remove all other equal to ia 
            if(print_bond(ia)) THEN 
               do ib = ia , nbugs
                  iret=bugs(ib) 
                  b1=ltor(1,iret)
                  b2=ltor(2,iret)
                  b3=ltor(3,iret)
                  b4=ltor(4,iret)
                  types2 = betb(b1)//betb(b2)//betb(b3)//betb(b4)
                  if(types1.eq.types2) print_bond(ib) = .false.
               end do
               b1=mend(nres(k,1)) 
               b2=mend(nres(l,1)) 
               b3=mend(nres(m,1)) 
               b4=mend(nres(n,1)) 
               write(kprint,1011) k,m,l,n, betb(k)
     &              ,betb(m),betb(l),betb(n),beta(k),beta(m),beta(l)
     &              ,beta(n),nres(k,1),prsymb(b1),nres(m,1),prsymb(b2)
     &              ,nres(l,1),prsymb(b3),nres(n,1),prsymb(b4)
1011           format(" *PTO",4i7,2x,' |',4(1x,a4),'|',4(1x,a4),'|',4(i4
     &              ,a3))      
            END IF
         END DO 
         IF(nbugs.eq.1001) THEN 
            errmsg = ' ..... and even more!!' 
            call xerror(errmsg,80,1,222)
         END IF
         call xerror(errmsg,80,1,223)
         notfound=1
      END IF

*=======================================================================
*---- Pick improper torsions that match input --------------------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Pick I-Torsions that Match Input   ---->'
      CALL pckitor(nbun,mend,natop,beta,jitor,nitor,NOATO,litr,litor,
     x            iret,errmsg)
      IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*---- Add extra Improprer Torsions -------------------------------------
*=======================================================================

      IF(mesg .EQ. 'SOLUTE:  ') THEN
         WRITE(kprint,'(5x,2a)')
     &        mesg,'Add extra I-Torsions               ---->'
         CALL addtor(beta,atres,ritor,xitor,xnitor,litr,litor,MM4,
     x        .TRUE.,concta,MM1,iret,errmsg)
         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*-------------- Match the torsion potential parameters -----------------
*----------- Iret is equal to 1 if the match is incomplete -------------
*=======================================================================

      WRITE(kprint,'(5x,2a)')
     &     mesg,'Match I-Torsions with Parameters   ---->'
      CALL pitors(litr,litor,MM4,betb,alito,pito,lpito,ITORP,
     x           potit,3,bugs,nbugs)

      IF(debug_it) THEN
        WRITE(kprint,'(a)'
     &       )'======= IMPROPER TORSION TABLE ======= '

          DO 171 m=1,litor
              ntph=DINT(potit(m,2)+0.5d0)
              WRITE(kprint
     &             ,'(1h ,4(a4,2x),2x,4(i4,1x),f7.2,i5,5x,i4)'
     &             )betb(litr(1,m)),betb(litr(2,m)),betb(litr(3,m))
     &             ,betb(litr(4,m)),litr(1,m),litr(2,m),litr(3,m),litr(4
     &             ,m),potit(m,1),ntph,m
171       CONTINUE
      END IF

      if(nbugs.ne.0)  THEN 
         call xerror(errmsg,80,1,221)
         errmsg = 'Improper torsion parameters not found for the'// 
     &        ' following quadruplets  - '
         call xerror(errmsg,80,1,222)
         errmsg = '               quadruplet        |'
     &        //'       types        |     labels         |  res     '
         call xerror(errmsg,80,1,222)
         do ia = 1,nbugs
            print_bond(ia)=.true.
         end do
c====    Start printing unassigned bonds 
         DO ia = 1,nbugs
            iret=bugs(ia) 
            k=litr(1,iret)
            m=litr(2,iret)
            l=litr(3,iret)
            n=litr(4,iret)
            types1 = betb(k)//betb(m)//betb(l)//betb(n)
c---        print i-tors ia and remove all lock print for all other 
            if(print_bond(ia)) THEN 
               do ib = ia , nbugs
                  iret=bugs(ib) 
                  b1=litr(1,iret)
                  b2=litr(2,iret)
                  b3=litr(3,iret)
                  b4=litr(4,iret)
                  types2 = betb(b1)//betb(b2)//betb(b3)//betb(b4)
                  if(types1.eq.types2) print_bond(ib) = .false.
               end do
               b1=mend(nres(k,1)) 
               b2=mend(nres(l,1)) 
               b3=mend(nres(m,1)) 
               b4=mend(nres(n,1)) 
               write(kprint,1012) k,m,l,n, betb(k)
     &              ,betb(m),betb(l),betb(n),beta(k),beta(m),beta(l)
     &              ,beta(n),nres(k,1),prsymb(b1),nres(m,1),prsymb(b2)
     &              ,nres(l,1),prsymb(b3),nres(n,1),prsymb(b4)
1012           format(" *PIT",4i7,2x,' |',4(1x,a4),'|',4(1x,a4),'|',4(i4
     &              ,a3))      
            END IF
         END DO 
         IF(nbugs.eq.1001) THEN 
            errmsg = ' ..... and even more!!' 
            call xerror(errmsg,80,1,222)
         END IF
         call xerror(errmsg,80,1,223)
         notfound=1
      END IF

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(notfound.eq.1) STOP

      IF(debug_it) THEN

*---------------- Print nhtype ----------------------------------------

        WRITE(kprint,*)
        WRITE(kprint,'(a)' )
     &       '=== HYDROGEN BOND INTEGER TABLE ==='
         DO 180 i=1,lpnbd
            DO 190 j=1,lpnbd
               IF(nhtype(i,j).NE.0) THEN
                  k=nhtype(i,j)
                  WRITE(kprint
     &                 ,'(2(a4,2x),i4,4x,2(i4,1x))' )alhyd(1,k),alhyd(2
     &                 ,k),k,i,j
               END IF
 190        CONTINUE
 180     CONTINUE
         
*--------------------------------------------------------------------
*----------- Print all the acceptor and donors of the protein ----------
*--------------------------------------------------------------------

         WRITE(kprint,*)
         WRITE(kprint,'(a)')
     &        '========== ACCEPTORS TABLE ========='
         DO 210 i=1,llacc
           WRITE(kprint,'(2(a4,2x),2(i4,2x),i4)' )
     x           beta(lacc(1,i)),beta(lacc(2,i)),lacc(1,i),lacc(2,i),
     x           i
 210     CONTINUE
         
         
         WRITE(kprint,*)
         WRITE(kprint,'(a)')
     &        '============ DONORS TABLE =========='
         DO 220 i=1,lldon
           WRITE(kprint,'(2(i4,2x),2(a4,2x),i4)' )
     x           ldon(1,i),ldon(2,i),beta(ldon(1,i)),beta(ldon(2,i) ),
     x           i
 220     CONTINUE
      END IF
      
      lconstr=0
      WRITE(kprint,10000) mesg,ntap,lbond,lbond,lconstr
     &     ,lbend,ltors,litor,int14p

*====================== END EXECUTABLE STATEMENTS ======================

10000 FORMAT(/
     &  '  *****************************************************',
     &  '*******************'/
     &  '  *                   ',a9,'  Initial TOPOLOGY  List',
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
20000 FORMAT(/
     &  '   =======================================================',
     &  '================='/
     &  '   =                                                      ',
     &  '                ='/
     &  '   =                 ',a9,' Assembling Molecules          ',
     &  '             ='/
     &  '   =                                                      ',
     &  '                ='/
     &  '   =======================================================',
     &  '================='/)
      RETURN
      END
