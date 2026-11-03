      SUBROUTINE rdatpg(kread,kprint,jerr,unit,nunit,n1,alpha,natype,
     x           charg1,natop,jngrp,jgrppt,cbond,nbond,cbend,nbend,
     x           ctor,ntor,citor,nitor,n3,term,back,n5,nback,jacc,
     x           jdon,nacc,ndon,n6,stpr,stsc,nstpr,nstsc,nstrg,n7,n8,
     x           iret,errmsg,mass1,ffchar)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine reads in the topology file.                      *
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
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  09/04/89 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS   Free format package, RDUNIT, CSORT, ATOMS,           *
*                 BOND, BACKBN, HBOND, NEAR0.                          *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n1,n3,n5,n6,n7,n8,iret,kread,kprint,jerr
      REAL*8 charg1(n3,*),mass1(n3,*)
      INTEGER nunit,nbun
      INTEGER natop(*),nbond(*),ntor(*),nitor(*),nback(*),nacc(*),
     x        ndon(*),jngrp(*),jgrppt(2,n3,*),nbend(*),nstpr(n7,*),
     x        nstsc(n7,*),nstrg(*)
      CHARACTER*7 alpha(n3,*),natype(n3,*),cbond(2,n3,*),ctor(4,n3,*),
     x            citor(4,n3,*),term(3,*),back(n5,*),jacc(2,n6,*),
     x            jdon(2,n6,*),cbend(3,n3,*),stpr(4,n7,*),stsc(n8,n7,*)
      CHARACTER*8 unit(*)
      CHARACTER*80 errmsg
      CHARACTER*1  ffchar
*-------------------- FUNCTION DEFINITION ------------------------------

      EXTERNAL near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,n,n2,m,nato,nex,nexa,nexd,nmc
      INTEGER ncter,index,ntape,nkeys,lna,ncount,ntmp
      INTEGER newkyv,newkyb,newkyc,newkyd,newkyf
      REAL*8  c12,c6
      LOGICAL error
      INTEGER nword,fioi,ilink,nprot(2,100),iline
      CHARACTER*8 alphh,char,space
      CHARACTER*80 line,strngs(40),err_atm
      CHARACTER*8 fmt,subcmd
      CHARACTER*26 err_str 
      CHARACTER*20 err_mis
      CHARACTER*23 err_unr
      CHARACTER*17 err_dim 
      CHARACTER*37 err_args(2)
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1,fior
      DATA sep/' ',','/comm/'#','!'/
      DATA space/'        '/

*==================== EXECUTABLE STATEMENTS ============================

      err_str = 'Error In TPG file at line '
      err_mis = ': Misplaced keyword '
      err_unr = ': Unrecognized keyword '
      err_dim = ' out of boundary.'
      err_args(1) = 'Number of arguments must not exceed  '
      err_args(2) = 'Number of arguments must be only     '
      line(79:80)='  '
      index=nunit
      iline = 0
200   READ(kread,'(a78)',END=600) line(1:78)
      iline = iline + 1
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 200 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         go to 200
      END IF

c=======================================================================
c              Start reading Residue 
c=======================================================================

      IF(strngs(1).EQ. 'RESIDUE' .OR. strngs(1) .EQ.
     x     'ENDGROUP' ) THEN
         char=strngs(2)(1:8)
         CALL usort(unit,n1,char,index,iret,errmsg)
         IF(iret.EQ.1) RETURN
         nbend(index)=0
         nback(index)=0
         nacc(index)=0
         ndon(index)=0
         nstrg(index)=0

 300     READ(kread,'(a78)',END=700) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 300 
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)

c=========keyword atoms ================================================

         IF(strngs(1).EQ. 'atoms' ) THEN
            j=0
            CALL atoms(kprint,kread,jerr,alpha(1,index),
     x           natype(1,index),charg1(1,index),jgrppt(1,1,index),
     x           jngrp(index),n3,nato,iline,j,errmsg,mass1(1,index))
            if(j.eq.1) THEN 
               iret=1
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg =err_str// fmt(1:j)
               RETURN
            END IF
            natop(index)=nato
            
c=========keyword rigid ================================================

         ELSE IF(strngs(1).EQ. 'rigid' ) THEN
            IF(natop(index) .EQ. 0) THEN
               j=0
               call int_str (iline,fmt,j)
               iret=1
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,3) 
               errmsg = err_str// fmt(1:j) //err_mis // '"rigid"'
               RETURN 
            ELSE
               nstrg(index)=nstrg(index)+1
               ncount=nstrg(index)
               
               IF(ncount .GT. n7) THEN
                  j=0
                  call int_str (iline,fmt,j)
                  iret=1
                  errmsg = char(1:8) // line(1:72) 
                  call xerror (errmsg,80,1,30) 
                  errmsg=
     &            err_str// fmt(1:j) //
     &                 ': Number of rigid units'//err_dim
                  RETURN
               END IF
               
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=24) ntmp
               CALL rdrgd(kprint,kread,jerr,ntmp,stpr(1,ncount,index),
     x              stsc(1,ncount,index),nstpr(ncount,index),
     x              nstsc(ncount,index),iret,iline,errmsg)
               IF(iret.EQ.1) THEN 
                  call xerror(errmsg,80,1,30) 
                  j=0
                  call int_str (iline,fmt,j)
                  errmsg=err_str// fmt(1:j) 
                  RETURN
               END IF

               IF(nstsc(ncount,index) .GT. n8) THEN
                  errmsg = char(1:8) // line(1:72) 
                  call xerror (errmsg,80,1,30) 
                  j=0
                  call int_str (iline,fmt,j)
                  errmsg=err_str//fmt(1:j)//': Secondary sites'//err_dim
                  iret=1
                  RETURN
               END IF

               IF(nstpr(ncount,index) .GT. 4 ) THEN
                  j=0
                  call int_str (iline,fmt,j)
                  errmsg=err_str// fmt(1:j)//': Primary sites'//err_dim
                  iret=1
                  RETURN
               END IF
            END IF

c=========keyword bonds ================================================

         ELSE IF(strngs(1).EQ. 'bonds' ) THEN
            subcmd='bonds'
            IF(natop(index) .EQ. 0) THEN
               j=0
               call int_str (iline,fmt,j)
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               errmsg=err_str// fmt(1:j)//err_mis// '"bonds"'
               iret=1
               RETURN
            END IF
            CALL bond(kprint,kread,jerr,cbond(1,1,index),n3,
     x           nmc,iret,iline,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF

            CALL chlab(subcmd,char,alpha(1,index),cbond(1,1,index),
     x           natop(index),nmc,2,iret,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            nbond(index)=nmc

c=========keyword omit_angles ==========================================

         ELSE IF(strngs(1).EQ. 'omit_angles' ) THEN
            subcmd='omit_angles'
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j)//err_mis//'"omit_angles"'
               iret=1
               RETURN
            END IF
            CALL bend(kprint,kread,jerr,cbend(1,1,index),n3,
     x           nmc,iret,iline,errmsg)
            IF(iret.EQ.1) THEN 
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF

            CALL chlab(subcmd,char,alpha(1,index),cbend(1,1,index),
     x           natop(index),nmc,3,iret,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            nbend(index)=nmc

c=========keyword dihed=================================================

         ELSE IF(strngs(1).EQ. 'dihed' ) THEN
            subcmd='dihed'
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j)//err_mis//'"dihed"'
               RETURN
            END IF
            CALL tors(kprint,kread,jerr,ctor(1,1,index),n3,
     x           nmc,iret,iline,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            CALL chlab(subcmd,char,alpha(1,index),ctor(1,1,index),
     x           natop(index),nmc,4,iret,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            ntor(index)=nmc

c=========keyword imphd=================================================

         ELSE IF(strngs(1).EQ. 'imphd' ) THEN
            subcmd='imphd'
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str//fmt(1:j)//err_mis // '"imphd"'
               RETURN
            END IF
            CALL tors(kprint,kread,jerr,citor(1,1,index),n3,
     x           nmc,iret,iline,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            CALL chlab(subcmd,char,alpha(1,index),citor(1,1,index),
     x           natop(index),nmc,4,iret,errmsg)
            IF(iret.EQ.1) THEN 
               call xerror(errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg=err_str// fmt(1:j) 
               RETURN
            END IF
            nitor(index)=nmc

c=========keyword backbone==============================================

         ELSE IF(strngs(1).EQ. 'backbone' ) THEN
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg =err_str// fmt(1:j)//err_mis//'"backbone"'
               RETURN
            END IF
            DO 40 i=2,nword
               back(nback(index)+i-1,index)=strngs(i)(1:7)
 40         CONTINUE
            nback(index)=nback(index)+nword-1

c=========keyword termatom==============================================

         ELSE IF(strngs(1).EQ. 'termatom' ) THEN
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg = err_str// fmt(1:j) // err_mis // '"termatom"'
            END IF
            term(1,index)(1:7)=strngs(2)(1:7)
            term(2,index)(1:7)=strngs(3)(1:7)
            IF(nword.EQ.4) THEN
               term(3,index)(1:7)=strngs(4)(1:7)
            ELSE
               term(3,index)=' '
            END IF

c=========keyword acc===================================================

         ELSE IF(strngs(1).EQ. 'acc' ) THEN
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg = err_str// fmt(1:j) // err_mis // '"acc"'
               RETURN
            END IF
            IF(nword .EQ. 3) THEN
               jacc(1,nacc(index)+1,index)=strngs(2)(1:7)
               jacc(2,nacc(index)+1,index)=strngs(3)(1:7)
               nacc(index)=nacc(index)+1
            ELSE IF(nword .EQ. 2) THEN
               jacc(1,nacc(index)+1,index)=strngs(2)(1:7)
               jacc(2,nacc(index)+1,index)=strngs(2)(1:7)
               nacc(index)=nacc(index)+1
            ELSE
               iret=1
               errmsg=err_args(1)//'1'
               RETURN
            END IF
               

c=========keyword don ==================================================

         ELSE IF(strngs(1).EQ. 'don' ) THEN
            IF(natop(index) .EQ. 0) THEN
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               iret=1
               j=0
               call int_str (iline,fmt,j)
               errmsg = err_str// fmt(1:j) // err_mis // '"don"'
               RETURN
            END IF
            IF(nword .EQ. 3) THEN
               jdon(1,ndon(index)+1,index)=strngs(2)(1:7)
               jdon(2,ndon(index)+1,index)=strngs(3)(1:7)
               ndon(index)=ndon(index)+1
            ELSE
               iret=1
               errmsg=err_args(2)//'1'
               RETURN
            END IF
               

c=========keyword blank line============================================

         ELSE IF(strngs(1) .EQ. ' ') THEN
            CONTINUE

c=========keyword RESIDUE_END ==========================================

         ELSE IF(strngs(1) .EQ. 'RESIDUE_END' .OR. strngs(1) .EQ. 'END'
     &           ) THEN

c------       Check termatoms

            subcmd='termatom'
            CALL chlab(subcmd,char,alpha(1,index),term(1,index),
     x           natop(index),1,2,iret,errmsg)
            IF(iret.eq.1) THEN 
               call xerror(errmsg,80,1,30)
               iret=0
            ENDIF   

c------       Check acceptors

            subcmd='acc'
            CALL chlab(subcmd,char,alpha(1,index),jacc(1,1,index),
     x           natop(index),nacc(index),2,iret,errmsg)
            IF(iret.eq.1) THEN 
               call xerror(errmsg,80,1,30)
               iret=0
            ENDIF   

c------       Check donors

            subcmd='don'
            CALL chlab(subcmd,char,alpha(1,index),jdon(1,1,index),
     x           natop(index),ndon(index),2,iret,errmsg)
            IF(iret.eq.1) THEN 
               call xerror(errmsg,80,1,30)
            ENDIF   

            IF(iret.EQ.1) THEN 
               errmsg = char(1:8) // line(1:72) 
               call xerror (errmsg,80,1,30) 
               j=0
               call int_str (iline,fmt,j)
               errmsg =err_str// fmt(1:j) 
               RETURN
            END IF
            
            GOTO 200

c=========Unknown Keyword ============================================

         ELSE
            errmsg = char(1:8) // line(1:72) 
            call xerror (errmsg,80,1,30) 
            j=0
            call int_str (iline,fmt,j)
            errmsg = err_str// fmt(1:j)//err_unr//strngs(1)(1:15)
            GOTO 500
         END IF
         GOTO 300
      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
      ELSE IF(strngs(1).EQ. 'END') THEN
         GOTO 600
         
      ELSE
         errmsg = char(1:8) // line(1:72) 
         call xerror (errmsg,80,1,30) 
         j=0
         call int_str (iline,fmt,j)
         if(strngs(1)(1:1).EQ.'&'.and.jerr.eq.1) THEN 
            errmsg = 'Missing termination keyword for residue '
     &           //char(1:8)  
         ELSE
            errmsg = err_str// fmt(1:j) //err_unr // strngs(1)(1:15)
         END IF
         GOTO 500
      END IF
      GOTO 200
      
 600  CONTINUE
      ! prepend FF char to second types from ligand file if required 
      if(ffchar.eq."T") THEN 
        do i=nunit+1,index
          do j=1,natop(i) 
            natype(j,i)(1:7)="*"//natype(j,i)(1:6)
          end do
        end do
      end if
      nunit=index
      iret=0
      RETURN
      
 700  iret=1
      j=0
      call int_str (iline,fmt,j)
      errmsg =err_str// fmt(1:j)//': Unexpected End Of File'
      RETURN
      
 500  iret=1
      RETURN 

 24   CONTINUE
      errmsg = char(1:8) // line(1:72) 
      call xerror (errmsg,80,1,30) 
      iret=1
      j=0
      call int_str (iline,fmt,j)
      errmsg = err_str // fmt(1:j) //': Internal reading error'
      call xerror (errmsg,80,1,30) 
      RETURN

      END


