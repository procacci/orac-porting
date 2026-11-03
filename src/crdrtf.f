      SUBROUTINE crdrtf(alpha,natype,charg1,natop,jngrp,jgrppt,cbond
     &     ,nbond,cbend,nbend,ctor,ntor,citor,nitor,k3,term,back,k5
     &     ,nback,jacc,jdon,nacc,ndon,k6,albon,pbon1,pbon2,lpbon,albnd
     &     ,pbnd1,pbnd2,pbnd3,pbnd4,pbnd5,pbnd6,pbnd7,lpbnd,altor,ptor1
     &     ,ntor2,potj,k7,lptor,alito,pito1,pito2,pito3,lpito,massu
     &     ,alhyd,xbond,xnbond,rbond,xtor,xntor,rtor,xitor,xnitor,ritor
     &     ,stpr,stsc,nstpr,nstsc,nstrg,k8,k9,iret,errmsg,mass1)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine reads in the topology file and the               *
*     potential parameter list for the run. There are 12               *
*     physical dimensions of the arrays in argument to be              *
*     provided. The subroutine contains two common statements          *
*     shared with the free format input package.                       *
*                                                                      *
*                                                                      *
*================== TOPOLOGY OF THE MACROMOLECULE =====================*
*                                                                      *
*     UNIT    :  List of labels for all the residues.                  *
*                >> character*8 UNIT(N1) <<                            *
*     NUNIT   :  Number of labels in the list of residues.             *
*     N1      :  Physical dimension of UNIT.                           *
*     ALPHA   :  List of all the atoms forming the residue.            *
*                >> character*7 ALPHA(N3,N2) <<                        *
*     NATYPE  :  List of the atom types.                               *
*                >> character*7 NATYPE(N3,N2) <<                       *
*     CHARG1  :  List of atomic charges. These charges are             *
*                used in the computation of the non-bonded             *
*                intramolecular interactions.                          *
*                >> real*8  CHARG1(N3,N2) <<                           *
*     NATOP   :  Number of atoms for all the residues.                 *
*                >> integer NATOP(N2) <<                               *
*     N2      :  Physical column dimension of the above                *
*                arrays.                                               *
*     CBOND    :  List of bonds for all the residues.                  *
*                >> character*7 CBOND(N3,N2) <<                        *
*     NBOND   :  Number of bond in the bond list.                      *
*                >> integer NBOND(N2) <<                               *
*     N3      :  Physical row dimension of the above arrays.           *
*     CORR    :  List of solvent-protein correlations.                 *
*                >> character*7 CORR(N4,N2) <<                         *
*     NCORR   :  Number of correlations.                               *
*                >> integer NCORR(N2) <<                               *
*     N4      :  Physical row dimension of CORR.                       *
*     TERM    :  List of terminal atoms. Maximum 3.                    *
*                >> character*7 TERM(3,N2) <<                          *
*     BACK    :  List of backbone atoms.                               *
*                >> character*7 BACK(N5,N2) <<                         *
*     N5      :  Physical row dimension of BACK.                       *
*     NBACK   :  Number of backbone atoms.                             *
*                >> integer NBACK(N2) <<                               *
*     JACC    :  List of acceptors.                                    *
*                >> character*7 JACC(N6,N2) <<                         *
*     JDON    :  List of donors.                                       *
*                >> character*7 JDON(N6,N2) <<                         *
*     NACC    :  Number of acceptors.                                  *
*                >> integer NACC(N2) <<                                *
*     NDON    :  Number of donors.                                     *
*                >> integer NDON(N2) <<                                *
*     N6      :  Physical row dimension of JACC/JDON.                  *
*                                                                      *
*===================== POTENTIAL PARAMETERS ===========================*
*                                                                      *
*     ALBND   :  List of bendings in the model potential.              *
*                >> character*7 ALBND(3,M1) <<                         *
*     PBND1   :  List of bending potential parameters                  *
*     PBND2      >> real*8 PBND1(M1),PBND2(M1) <<                      *
*     LPBND   :  Number of bendings in the model potential.            *
*     M1      :  Physical column dimension of ALBND.                   *
*     ALTOR   :  List of proper torsions in the model potential.       *
*                >> character*7 ALTOR(4,M2) <<                         *
*     PTOR1   :  List of potential parameters for proper               *
*     NTOR2      torsions.                                             *
*                >> real*8 PTOR1(M2), integer NTOR2(M2) <<             *
*     LPTOR   :  Number of proper torsions in the model                *
*                potential.                                            *
*     M2      :  Physical column dimension of ALTOR.                   *
*     ALITO   :  List of improper torsions in the model potential.     *
*                >> character*7 ALITO(4,M3) <<                         *
*     PITO1   :  List of potential parameters for improper             *
*     PITO2      torsions.                                             *
*                >> real*8 PITO1(M3), integer PITO2(M3) <<             *
*     LPITO   :  Number of improper torsions in the model              *
*                potential.                                            *
*     M3      :  Physical column dimension of ALITO.                   *
*     ALNBD   :  List of non-bonded interactions in the                *
*                model potential.                                      *
*                >> character*7 ALNBD(*) <<                            *
*     PNBD1   :  List of non-bonded potential parameters.              *
*     PNBD2      >> real*8 PNBD1(*), PNBD2(*), PNBD3(*) <<             *
*     PNBD3      Input in Kcal/mol.                                    *
*     MASS    :  List of atomic masses.                                *
*             :  >> real*8 MASS(*) <<                                  *
*     LPNBD   :  Number of non-bonded interactions in the              *
*                model potential.                                      *
*     ECC12   :  List of non-bonded potential parameters               *
*     ECC6       read when the mixing rules are not used.              *
*                >> real*8 ECC12(*), ECC6(*) <<                        *
*                In Kcal/mol.                                          *
*     IZ      :  Mixing rules flag. IZ=0 mixing rules are used.        *
*                IZ=1 mixing rules are not used.                       *
*     ALHYD   :  List of hydrogen bonds in the model                   *
*                potential.                                            *
*                >> character*7 ALHYD(2,*) <<                          *
*     PHYD1   :  List of hydrogen bond potential parameters.           *
*     PHYD2      >> real*8  PHYD1(M5), PHYD2(M5) <<                    *
*     LPHYD   :  Number of hydrogen bond parameters in the             *
*                model potential.                                      *
*     M5      :  Physical dimension of PHYD1/PHYD2.                    *
*                                                                      *
*=================== LIST OF RESIDUES FOR THE PROTEIN =================*
*                                                                      *
*     JOIN    :  List of residues occurring in the actual              *
*                protein.                                              *
*                >> integer JOIN(M6) <<                                *
*     M6      :  Physical dimension of JOIN                            *
*     NBUN    :  Number of residues constituting the actual            *
*                protein.                                              *
*                                                                      *
*=================== LIST OF EXTRA LINKS ==============================*
*                                                                      *
*     LINK    :  List of extra links for the solute atoms.             *
*                >> integer LINK(2,*) <<                               *
*     NLINK   :  Total number of extra links.                          *
*                                                                      *
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*================== Description of the keywords =======================*
*                                                                      *
* READ_PFR_BIN            filename                                     *
* WRITE_PFR_BIN           filename                                     *
* READ_TPF_ASCII          filename {PRINT}                             *
* READ_PRM_ASCII          filename {PRINT}                             *
* REPL_RESIDUE            $                                            *
*                         (Accept same keywords as RESIDUE environment)*
* REPL_RESIDUE_END        $                                            *
* PRINT_TOPOLOGY          string [a;b;c;i;h;n;4;p;w]                   *
* JOIN                    $                                            *
*                 frag_begin      $                                    *
*                         [char;......]  "{4 x char = char char        *
*                                    char char}"                       *
*                 frag_end $                                           *
* JOIN_END                                                             *
* ADD_TPG                 $                                            *
*                         (Accept same keywords as RESIDUE environment)*
* ADD_TPG_END             $                                            *
* CHROMOPHORE             $                                            *
*                 chromo          $                                    *
*                         group                                        *
*                                 .....                                *
*                                 ffield(nfatm(nfield)++,nfield) ....  *
*                                 .....                                *
*         (can use format n1 - n2 to take all atoms between n1 and n2) *
*                         end $                                        *
*                         exclude                                      *
*                                 .....                                *
*                                 ffield(nfatm(nfield)++,nfield) ....  *
*                                 .....                                *
*         (can use format n1 - n2 to take all atoms between n1 and n2) *
*                         end $                                        *
*                 end             $                                    *
* CHROMOPHORE_END         $                                            *
* READ_FIELD              mfield  filename                             *
* START_RUN                                                            *
*                                                                      *
*================== Description of the keywords end ===================*
*                                                                      *
*     $Id: crdrtf.f,v 1.1.1.1 2009/10/19 07:48:42 signo Exp $
*                                                                      *
*     Written by Massimo Marchi CE Saclay 1994                         *
*     Updated by P Procacci at CECAM Jan 1997                          *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER k3,k5,k6,k7,k8,k9,iret
      REAL*8 massu(*),charg1(k3,*),mass1(k3,*),pbnd1(*),pbnd2(*),pbnd3(
     &     *),pbnd4(*),pbnd5(*),pbnd6(*),pbnd7(*),ptor1(*),pito1(*)
     &     ,pito2(*),pito3(*),potj(k7,*),pbon1(*),pbon2(*)
      INTEGER lpbnd,lptor,lpito,xnbond,xntor,xnitor,lpbon
      INTEGER natop(*),nbond(*),ntor(*),nitor(*),nback(*),nacc(*),
     x        ndon(*),jngrp(*),ntor2(*),nbend(*),rbond(2,*),rtor(4,*),
     x        ritor(4,*),jgrppt(2,k3,*),nstpr(k8,*),nstsc(k8,*),
     x        nstrg(*)
      CHARACTER*7 alpha(k3,*),natype(k3,*),cbond(2,k3,*),ctor(4,k3,*),
     x            citor(4,k3,*),term(3,*),back(k5,*),jacc(2,k6,*),
     x            jdon(2,k6,*),xbond(2,*),xtor(4,*),xitor(4,*),
     x            cbend(3,k3,*),stpr(4,k8,*),stsc(k9,k8,*)
      CHARACTER*7 albnd(3,*),altor(4,*),alito(4,*),alhyd(2,*),albon(2,*)
      CHARACTER*80 errmsg,message0,message1

*-------------------- FUNCTION DEFINITION ------------------------------
      EXTERNAL near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,j,n,m,o1,kdummy,nsevere,nwarning
      INTEGER index,ntape,totdf
      INTEGER flag
      LOGICAL bin0,bin1,f_ex,ok,chrflg,warn_bin,join_solvent
     &     ,join_solute,ok1,ok2
      INTEGER nword,ilink
      CHARACTER*8 char,space
      CHARACTER*80 line,strngs(40),chartmp
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2),ffchar
      DATA sep/' ',','/comm/'#','!'/


*-------------------- COMMON VARIABLES ---------------------------------


      COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
      INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)
      DATA space/'        '/

*==================== EXECUTABLE STATEMENTS ============================

      nsevere= 0
      nwarning = 0
      ilink=0
      nunit=0
      nbun=0
      nbun_slv=0
      bin0=.FALSE.
      bin1=.FALSE.
      chrflg=.FALSE.
      join_solute=.FALSE.
      join_solvent=.FALSE.
      ok1=.FALSE.
      ok2=.FALSE.
      lpbon=0    
      lpbnd=0
      lptor=0
      lpito=0
      lpnbd=0
      lphyd=0
      xnbond=0
      xnitor=0
      line(79:80)='  '
      warn_bin=.false.
      ffchar="1" 

c=======================================================================
c     Environment parser starts here: This is &PARAMETERS env. 
c=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Command  READ_TPGPRM_BIN=========================================
c---- Read potential parameters and topology from an unformatted ------
c---- file.

      IF(strngs(1).EQ. 'READ_TPGPRM_BIN'.OR.strngs(1).EQ. 'READ_PFR_BIN'
     &     ) THEN
         bin0=.TRUE.
         tpgfil = .TRUE.
         IF(.NOT. tpgwbn) THEN
            IF(nword .LT. 2) THEN
               iret=1
               errmsg=' Error in CRDRTF: Provide directive '//
     x              'READ_TPGPRM_BIN with Filename. Abort. '
               RETURN
            END IF
            CALL uscrpl(strngs(2),80)
            INQUIRE(FILE=strngs(2),EXIST=f_ex)
            IF(f_ex) THEN
               CALL openf(ktpgprm_read,strngs(2),'FORMATTED','OLD',0)
            ELSE
               errmsg=' Error in CRDRTF: File does not exist. Abort '
               iret=1
               RETURN
            END IF
         END IF
         
c=====Command WRITE_TPGPRM_BIN ============================================
c----- Write potential parameters and topology to an unformatted -------
c----- file. -----------------------------------------------------------

       ELSE IF(strngs(1).EQ. 'WRITE_TPGPRM_BIN'.OR.strngs(1).EQ
     &        .'WRITE_PFR_BIN') THEN
         tpgwbn = .TRUE.
         IF(nword .LT. 2) THEN
            nsevere = nsevere + 1
            errmsg=' Missing file name' 
            CALL xerror(errmsg,80,1,30)
         ELSE
            CALL uscrpl(strngs(2),80)
            INQUIRE(FILE=strngs(2),EXIST=f_ex)
            IF(f_ex) THEN
               CALL openf(ktpgprm_write,strngs(2),'FORMATTED','OLD'
     &              ,0)
            ELSE
               CALL openf(ktpgprm_write,strngs(2),'FORMATTED','NEW'
     &              ,0)
            END IF
         END IF

c=====Command READ_TPG_ASCII ===========================================
c-----Read topology from a text file ----------------------------------

       ELSE IF(strngs(1).EQ. 'READ_TPG_ASCII'.OR.strngs(1).EQ.
     &        'READ_TPF_ASCII') THEN
         IF(nword .LT. 2) THEN
            nsevere = nsevere + 1
            errmsg=' Missing file name: Provide topology file name.' 
            CALL xerror(errmsg,80,1,30)
         ELSE
c--         Topology file being read in ...            
            if(ok1) THEN
              ffchar="T"
            else
              ok1=.TRUE.
            end if
            if(join_solute) THEN
              nsevere = nsevere + 1
              errmsg=' All topology must be given BEFORE JOIN_SOLUTE' 
              CALL xerror(errmsg,80,1,30)
            ELSE
              bin1=.TRUE.
            END IF
            CALL uscrpl(strngs(2),80)
            INQUIRE(FILE=strngs(2),EXIST=f_ex)
            if(f_ex)  THEN 
               CALL openf(ntape,strngs(2),'FORMATTED','OLD',0)
               IF(strngs(3) .EQ. 'PRINT') THEN
                  flag=1
               ELSE
                  flag=0
               END IF
               CALL rdatpg(ntape,kprint,flag,prsymb,nunit,n5,alpha,
     x              natype,charg1,natop,jngrp,jgrppt,cbond,nbond,cbend,
     x              nbend,ctor,ntor,citor,nitor,k3,term,back,k5,nback,
     x              jacc,jdon,nacc,ndon,k6,stpr,stsc,nstpr,nstsc,
     x              nstrg,k8,k9,iret,errmsg,mass1,ffchar)
               IF(iret .EQ. 1) THEN 
                  call xerror(errmsg,80,1,30) 
                  errmsg = ' 1 ERROR WHILE EXECUTING RDATPG' 
                  call xerror(errmsg,80,1,2) 
                  STOP
                ENDIF 
            ELSE
               errmsg='Topology File not found. ABORT'
               CALL xerror(errmsg,80,1,2)
               STOP
            END IF
         END IF

c=====Command READ_PRM_ASCII============================================
c---- Read potential parameters from a text file ----------------------

      ELSE IF(strngs(1).EQ. 'READ_PRM_ASCII') THEN
         IF(nword .LT. 2) THEN
            nsevere = nsevere + 1
            errmsg=' Missing file name: Provide parameter file name.' 
            CALL xerror(errmsg,80,1,30)
         ELSE
            if(ok2) THEN
              ffchar="R"
            else
              ok2=.TRUE.
            end if
            CALL uscrpl(strngs(2),80)
            INQUIRE(FILE=strngs(2),EXIST=f_ex)
            if(f_ex)  THEN 
               CALL openf(ntape,strngs(2),'FORMATTED','OLD',0)
               IF(strngs(3) .EQ. 'PRINT') THEN
                  flag=1
               ELSE
                  flag=0
               END IF
               CALL rdaprm(ntape,kprint,flag,albon,pbon1,pbon2,lpbon
     &              ,albnd,pbnd1,pbnd2,pbnd3,pbnd4,pbnd5,pbnd6,pbnd7
     &              ,lpbnd,altor,ptor1,ntor2,potj,k7,lptor,alito,pito1
     &              ,pito2,pito3,lpito,alnbd,pnbd1,pnbd2,pnbd3,pnbd4,m6
     &              ,ejorg,sjorg,lpnbd,massu,ecc12,ecc6,iz,alhyd,phyd1
     &              ,phyd2,lphyd,ffchar,iret,errmsg)
               IF(iret .EQ. 1) THEN 
                  call xerror(errmsg,80,1,30) 
                  errmsg = ' ERROR WHILE EXECUTING RDAPRM' 
                  call xerror(errmsg,80,1,2) 
                  STOP
               END IF
            ELSE
               errmsg='Parameter File not found. ABORT'
               CALL xerror(errmsg,80,1,2)
               STOP
            END IF
         END IF
c=====Command REPL_RESIDUE==============================================
c---- Add or replace a new residue to the list ------------------------

      ELSE IF(strngs(1).EQ. 'REPL_RESIDUE') THEN
         IF (bin1) THEN 
            flag=0
            CALL rdatpg(knlist,kprint,flag,prsymb,nunit,n5,alpha,
     x           natype,charg1,natop,jngrp,jgrppt,cbond,nbond,cbend,
     x           nbend,ctor,ntor,citor,nitor,k3,term,back,k5,nback,
     x           jacc,jdon,nacc,ndon,k6,stpr,stsc,nstpr,nstsc,
     x           nstrg,k8,k9,iret,errmsg,mass1,ffchar)
            go to 100
            IF(iret .EQ. 1) THEN 
               call xerror(errmsg,80,1,30) 
               errmsg = ' ERROR WHILE EXECUTING RDATPG' 
               call xerror(errmsg,80,1,2) 
               STOP
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg='Topology unknown: specify READ_TPG_ASCII before' //
     &           ' REPL_RESIDUE' 
            CALL xerror(errmsg,80,1,30)
            GO TO 600 
         END IF

c=====Command PRINT_TOPOLOGY ===========================================
c----- Print topology and potential parameters for the run -------------

      ELSE IF(strngs(1).EQ. 'PRINT_TOPOLOGY') THEN
         prttpg=.TRUE.
900      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 900
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
      
         IF(strngs(1) .EQ. 'bonds') THEN
            prtbndl=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'sequence') THEN
            prtseq=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'constraints') THEN
            prtcnl=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'atoms') THEN
            prtatl=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'bendings') THEN
            prtbal=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'P-torsions') THEN
            prtptl=.TRUE.
            
         ELSE IF(strngs(1) .EQ. 'I-torsions') THEN
            prtitl=.TRUE.

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 900

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg='Unrecognized keyword ----> '
     &           //strngs(2)//'...or missing &END'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 900

c=====Command JOIN======================================================
c----- Join the protein residue together -------------------------------

      ELSE IF(strngs(1).EQ. 'JOIN' ) THEN
         IF(strngs(2) .EQ. ' ' .OR. strngs(2) .EQ. 'SOLUTE') THEN
            IF(join_solute) THEN
               nsevere = nsevere + 1
               errmsg='Can call JOIN SOLUTE only once.'
               CALL xerror(errmsg,80,1,30)
               GO TO 600
            END IF
            join_solute=.TRUE.
            slt_exist=.TRUE.
            IF (bin1) THEN
               m=0
               n=0
               nplmap=0
800            READ(knlist,'(a78)',END=600) line(1:78)
               IF(jerr.EQ.1) CALL wrenc(kprint,line)
               IF(line(1:1) .EQ. '#') GOTO 800 
               CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
               IF(nword.EQ.0) GOTO 800
               IF(strngs(1) .NE. 'END' ) THEN
                  ok=.TRUE.
                  DO i=1,nword
                     IF(.NOT.ok) GOTO 1850
                     char=strngs(i)(1:8)
                     IF(char .NE. 'x' .AND. char .NE. 'X') THEN
                        CALL csort(prsymb,nunit,char,index,iret,
     x                       errmsg)
                        if(iret.eq.1) THEN 
                           nsevere = nsevere + 1
                           CALL xerror(errmsg,80,1,30)
                           iret=0
                        ELSE
                           n=n+1
                           IF(n .GT. nores) THEN
                              i1 = nores
                              call int_str(i1,fmt,j)
                              errmsg= fmt(1:j) // ' residues defined: '/
     &                             /'physical dimensions exceed'
                              nsevere = nsevere + 1
                              CALL xerror(errmsg,80,1,30)
                              go to 600
                           ELSE
                              mend(n)=index
                           END IF
                        ENDIF
                     ELSE
                        char=strngs(i-1)(1:8)
                        CALL csort(prsymb,nunit,char,index
     &                       ,iret,errmsg)
                        if(iret.eq.1) THEN 
                           nsevere = nsevere + 1
                           CALL xerror(errmsg,80,1,30)
                           iret=0
                        ELSE
                           CALL fndfmt(1,strngs(i+1),fmt)
                           READ(strngs(i+1),fmt,err=20) o1
                           DO j=1,o1-1
                              n=n+1
                              IF(n .GT. nores) THEN
                                 i1 = nores
                                 call int_str(i1,fmt,j)
                                 errmsg=
     &                                fmt(1:j) // ' residues defined: '/
     &                                /'physical dimensions exceed'
                                 nsevere = nsevere + 1
                                 CALL xerror(errmsg,80,1,30)
                                 go to 600
                              ELSE
                                 mend(n)=index
                              END IF
                           END DO
                           ok=.FALSE.
                           GOTO 1860
                        END IF
                     END IF
1850                 CONTINUE
                     ok=.TRUE.
1860                 CONTINUE
                  END DO
                  GOTO 800
               END IF
               nbun=n
            ELSE IF(tpgfil) THEN 
200            READ(knlist,'(a78)',END=600) line(1:78)
               IF(jerr.EQ.1) CALL wrenc(kprint,line)
               IF(line(1:1) .EQ. '#') GOTO 200 
               CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x              errmsg)
               IF(iret.EQ.0) THEN
                  IF(strngs(1).NE. 'END') GOTO 200
                  GOTO 100
               END IF

            ELSE
               nsevere = nsevere + 1
               errmsg='Topology unknown: specify READ_TPG_ASCII before '
     &              //'JOIN' 
               CALL xerror(errmsg,80,1,30)
               GO TO 600
            END IF
         ELSE IF(strngs(2) .EQ. 'SOLVENT') THEN 
            IF(join_solvent) THEN
               nsevere = nsevere + 1
               errmsg='Can call JOIN SOLVENT only once.'
               CALL xerror(errmsg,80,1,30)
               GO TO 600
            END IF
            join_solvent=.TRUE.
            IF (bin1) THEN 
               m=0
               n=0
               nplmap=0
810            READ(knlist,'(a78)',END=600) line(1:78)
               IF(jerr.EQ.1) CALL wrenc(kprint,line)
               IF(line(1:1) .EQ. '#') GOTO 810 
               CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
               IF(nword.EQ.0) GOTO 810
               IF(strngs(1) .NE. 'END' ) THEN
                  ok=.TRUE.
                  DO i=1,nword
                     IF(.NOT.ok) GOTO 1851
                     char=strngs(i)(1:8)
                     IF(char .NE. 'x' .AND. char .NE. 'X') THEN
                        CALL csort(prsymb,nunit,char,index,iret,
     x                       errmsg)
                        if(iret.eq.1) THEN 
                           nsevere = nsevere + 1
                           CALL xerror(errmsg,80,1,30)
                           iret=0
                        ELSE
                           n=n+1
                           IF(n .GT. nores) THEN
                              i1 = nores
                              call int_str(i1,fmt,j)
                              errmsg= fmt(1:j) // ' residues defined: '/
     &                             /'physical dimensions exceed'
                              nsevere = nsevere + 1
                              CALL xerror(errmsg,80,1,30)
                              go to 600
                           ELSE
                              mend_slv(n)=index
                           END IF
                        ENDIF
                     ELSE
                        char=strngs(i-1)(1:8)
                        CALL csort(prsymb,nunit,char,index
     &                       ,iret,errmsg)
                        if(iret.eq.1) THEN 
                           nsevere = nsevere + 1
                           CALL xerror(errmsg,80,1,30)
                           iret=0
                        ELSE
                           CALL fndfmt(1,strngs(i+1),fmt)
                           READ(strngs(i+1),fmt,err=20) o1
                           DO j=1,o1-1
                              n=n+1
                              IF(n .GT. nores) THEN
                                 i1 = nores
                                 call int_str(i1,fmt,j)
                                 errmsg=
     &                                fmt(1:j) // ' residues defined: '/
     &                                /'physical dimensions exceed'
                                 nsevere = nsevere + 1
                                 CALL xerror(errmsg,80,1,30)
                                 go to 600
                              ELSE
                                 mend_slv(n)=index
                              END IF
                           END DO
                           ok=.FALSE.
                           GOTO 1861
                        END IF
                     END IF
1851                 CONTINUE
                     ok=.TRUE.
1861                 CONTINUE
                  END DO
                  GOTO 810
               END IF
               nbun_slv=n
            ELSE IF(tpgfil) THEN 
210            READ(knlist,'(a78)',END=600) line(1:78)
               IF(jerr.EQ.1) CALL wrenc(kprint,line)
               IF(line(1:1) .EQ. '#') GOTO 210
               CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x              errmsg)
               IF(iret.EQ.0) THEN
                  IF(strngs(1).NE. 'END') GOTO 210
                  GOTO 100
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg='Topology unknown: specify READ_TPG_ASCII before '
     &              //'JOIN' 
               CALL xerror(errmsg,80,1,30)
               GO TO 600
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg='Argument to JOIN incorrect: Must specify either '
     &           //' SOLUTE or SOLVENT' 
            CALL xerror(errmsg,80,1,30)
            GO TO 600
         END IF
      
c====Command ADD_TPG====================================================
c----- Add extra bonds, torsions and improper tosions ------------------
      
      ELSE IF(strngs(1) .EQ. 'ADD_TPG' .AND. strngs(2) .EQ. 'SOLUTE')
     &        THEN
         CALL addtpg(knlist,kprint,jerr,xbond,xnbond,rbond,
     x        xtor,xntor,rtor,xitor,xnitor,ritor,iret,errmsg)
         IF(iret .EQ. 1) THEN 
            call xerror(errmsg,80,1,30) 
            errmsg = ' ERROR WHILE EXECUTING ADDTPG' 
            call xerror(errmsg,80,1,2) 
            STOP
         END IF


c=====Command CHROMOPHORE =============================UNSUPPORTED======
c----- Define which region of the solute to consider as a chromophore---

      ELSE IF(strngs(1) .EQ. 'CHROMOPHORE') THEN
         errmsg = 'UNSUPPORTED COMMAND'   
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         CALL getchr(ffield,fexcl,nfatm,nexatm,nfield,f1,f2,f3,
     x        iret,errmsg)
         chrflg=.TRUE.
         IF(iret .EQ. 1) RETURN

c=====Command READ_FIELD================================UNSUPPORTED=====

      ELSE IF(strngs(1) .EQ. 'READ_FIELD') THEN
         errmsg = 'UNSUPPORTED COMMAND'   
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         IF(.NOT.chrflg) THEN
            nsevere = nsevere + 1
            errmsg='specify READ_FIELD after CHROMOPHORE '
            CALL xerror(errmsg,80,1,30)
         ELSE
            lfield=.TRUE.
            IF (nword .NE. 3) THEN
               nsevere = nsevere + 1
               errmsg='Number of arguments must be 2'
               CALL xerror(errmsg,80,1,30)
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) mfield
               totdf=0
               DO i=1,nfield
                  totdf=totdf+nfatm(i)
               END DO
               totdf=16*totdf+8
               chartmp='DATA$'//strngs(3)(1:75)
               CALL openf(kfield,strngs(3),'UNFORMATTED','UNKNOWN',0)
               CALL openf(kdummy,chartmp,'FORMATTED','UNKNOWN',0)
            END IF
         END IF

c====Missing &END Keyword ==============================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg='New ENV.' //strngs(1)(1:12)//
     &        ' begins before &END keyword.'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600


c=====&END Keyword : exit ==============================================

      ELSE IF(strngs(1).EQ. '&END' ) THEN
         GOTO 600

c=====Comand STOP_RUN ==================================================
c---- Do only input operations ----------------------------------------

      ELSE IF(strngs(1).EQ. 'STOP_RUN' ) THEN
         stoprun=.TRUE.

c=====Blank Line========================================================

      ELSE IF(strngs(1).EQ. ' ' ) THEN
         CONTINUE

c==== Wrong keyword ====================================================

      ELSE
         errmsg='unrecognized command '// strngs(1)//
     &        ' or missing &END'
         call xerror(errmsg,80,1,30)
         nsevere = nsevere + 1 
      END IF

      GOTO 100

600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

      IF(ok1 .AND. ok2) tpgprm_read=.TRUE.

      IF(warn_bin) THEN
         errmsg = 'Ignored Input Ended ' 
         CALL xerror(errmsg,80,1,11)
      END IF

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         call int_str(nsevere,fmt,j)
         goto 500 
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_PARAMETERS'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_PARAMETERS'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_PARAMETERS'
         call xerror(errmsg,80,1,1)
      ENDIF    

      RETURN

c==============================================================================
c     Errors were found
c==============================================================================

500   CONTINUE

      errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_PARAMETERS'
      CALL xerror(errmsg,80,1,2)
      RETURN

 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      RETURN

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
