      SUBROUTINE rdaprm(kread,kprint,jerr,albon,pbon1,pbon2,lpbon,albnd
     &     ,pbnd1,pbnd2,pbnd3,pbnd4,pbnd5,pbnd6,pbnd7,lpbnd,altor,ptor1
     &     ,ntor2,potj,k1,lptor,alito,pito1,pito2,pito3,lpito,alnbd
     &     ,pnbd1,pnbd2,pnbd3,pnbd4,atop,ejorg,sjorg,lpnbd,mass,ecc12
     &     ,ecc6,iz,alhyd,phyd1,phyd2,lphyd,ffchar,iret,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine reads in the potential parameters                *
*     list.                                                            *
*                                                                      *
*                                                                      *
*======================================================================*
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
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  03/01/92 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley CA 1992                    *
*                                                                      *
*     EXTERNALS   Free format package, RDUNIT, CSORT, ATOMS,           *
*                 BOND, BACKBN, HBOND, NEAR0.                          *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iz,iret,jerr,kread,kprint,k1
      REAL*8 pbnd1(*),pbnd2(*),pbnd3(*),pbnd4(*),pbnd5(*),pbnd6(*)
     &     ,pbnd7(*),ptor1(*),potj(k1,*),pito1(*),pito2(*),pito3(*)
     &     ,pnbd1(*),pnbd2(*),pnbd3(*),pnbd4(*),mass(*),ecc12(*),ecc6(*)
     &     ,phyd1(*),phyd2(*),ejorg(*),sjorg(*),pbon1(*),pbon2(*)
      INTEGER ntor2(*)
      INTEGER lpbnd,lptor,lpito,lpnbd,lphyd,lpbon,atop
      CHARACTER*7 albnd(3,*),altor(4,*),alito(4,*),alnbd(*),
     x            alhyd(2,*),albon(2,*)
      CHARACTER*80 errmsg
      CHARACTER*1  ffchar

*-------------------- FUNCTION DEFINITION ------------------------------

      EXTERNAL near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,j1,ij,iline,nsevere,k,l
      INTEGER nword
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      CHARACTER*26 err_str 
      CHARACTER*20 err_mis
      CHARACTER*23 err_unr
      CHARACTER*17 err_dim 
      CHARACTER*35 err_args
      REAL*8  v0,en,gamma,dummy
      INTEGER nmax_bond
      INTEGER nmax_bend
      INTEGER nmax_ptor
      INTEGER nmax_itor
      INTEGER loc_pbon,loc_pbnd,loc_ptor,loc_pito,loc_pnbd,loc_phyd

      COMMON /dim_for_prm/ nmax_bond,nmax_bend,nmax_ptor,nmax_itor
      DATA sep/' ',','/comm/'!','#'/

*==================== EXECUTABLE STATEMENTS ============================


      err_str = 'Error In PRM file at line '
      err_mis = ': Misplaced keyword '
      err_unr = ': Unrecognized keyword '
      err_dim = ' out of boundary.'
      err_args= ' parameters expected after keyword ' 
      line(79:80)='  '
      nsevere = 0
      iline = 0 

      loc_pbon=lpbon    
      loc_pbnd=lpbnd
      loc_ptor=lptor
      loc_pito=lpito
      if(iz.eq.0) loc_pnbd=lpnbd
      loc_phyd=loc_phyd


c=======================================================================
c     Parameters file parsing starts...
c=======================================================================


100   READ(kread,'(a78)',END=600) line(1:78)
      iline = iline + 1
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      iret=0 
      IF(iret.EQ.1) THEN 
         j=0
         k=iline
         call int_str (k,fmt,j)
         errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
         call xerror (errmsg,80,1,222) 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         go to 100
      END IF

c===  keyword BENDING ===================================================

      IF(strngs(1).EQ. 'BENDINGS'. OR. strngs(1).EQ. 'ANGLES' ) THEN
110      READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 110 
         iret=0 
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 110
         END IF
         IF(nword.EQ.0) GOTO 110
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 110
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF
         lpbnd=lpbnd+1

c---     Out of boundaries: STOP 
         IF(lpbnd .GT. nmax_bend) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='Bends'// err_dim//
     &           ' Change "bendp" in parameters.h'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF

c--      Error if args are not 5 
         IF(nword .EQ. 5) THEN
            albnd(1,lpbnd)(1:7)=strngs(1)(1:7)
            albnd(2,lpbnd)(1:7)=strngs(2)(1:7)
            albnd(3,lpbnd)(1:7)=strngs(3)(1:7)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pbnd1(lpbnd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pbnd2(lpbnd)
            pbnd3(lpbnd)=0.0D0
            pbnd4(lpbnd)=0.0D0
            GOTO 110
         ELSE IF(nword .EQ. 7) THEN            
            albnd(1,lpbnd)(1:7)=strngs(1)(1:7)
            albnd(2,lpbnd)(1:7)=strngs(2)(1:7)
            albnd(3,lpbnd)(1:7)=strngs(3)(1:7)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pbnd1(lpbnd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pbnd2(lpbnd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) pbnd3(lpbnd)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) pbnd4(lpbnd)
            GOTO 110
         ELSE IF(nword .EQ. 8) THEN            
            albnd(1,lpbnd)(1:7)=strngs(1)(1:7)
            albnd(2,lpbnd)(1:7)=strngs(2)(1:7)
            albnd(3,lpbnd)(1:7)=strngs(3)(1:7)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pbnd1(lpbnd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pbnd2(lpbnd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) pbnd3(lpbnd)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) pbnd4(lpbnd)
            CALL fndfmt(2,strngs(8),fmt)
            READ(strngs(8),fmt,err=20) pbnd5(lpbnd)
            GOTO 110
         ELSE IF(nword .EQ.10) THEN            
            albnd(1,lpbnd)(1:7)=strngs(1)(1:7)
            albnd(2,lpbnd)(1:7)=strngs(2)(1:7)
            albnd(3,lpbnd)(1:7)=strngs(3)(1:7)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pbnd1(lpbnd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pbnd2(lpbnd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) pbnd3(lpbnd)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) pbnd4(lpbnd)
            CALL fndfmt(2,strngs(8),fmt)
            READ(strngs(8),fmt,err=20) pbnd5(lpbnd)
            CALL fndfmt(2,strngs(9),fmt)
            READ(strngs(9),fmt,err=20) pbnd6(lpbnd)
            CALL fndfmt(2,strngs(10),fmt)
            READ(strngs(10),fmt,err=20) pbnd7(lpbnd)
            GOTO 110
         ELSE
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '5'// err_args // 'BENDING'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GOTO 110
         END IF

c=====keyword BOND  ==================================================
         
      ELSE IF(strngs(1).EQ. 'BOND'. OR. strngs(1).EQ. 'BONDS') THEN
170      READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#'.OR.line(1:1).EQ.'!') GOTO 170 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 170
         END IF
         IF(nword.EQ.0) GOTO 170
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 170
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF
         lpbon=lpbon+1

c---     Out of boundaries: STOP 
         IF(lpbon .GT. nmax_bond) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='Bonds'// err_dim//
     &           ' Change "bondp" in parameters.h'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF

         IF(nword .NE. 4) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '4'// err_args // 'BOND'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GO TO 170 
         ELSE
            albon(1,lpbon)(1:7)=strngs(1)(1:7)
            albon(2,lpbon)(1:7)=strngs(2)(1:7)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pbon1(lpbon)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pbon2(lpbon)
            GOTO 170
         END IF

c=====keyword TORSION ===============================================

      ELSE IF((strngs(1).EQ. 'TORSION' .AND. strngs(2) .EQ.
     x        'PROPER') .OR.  strngs(1).EQ. 'DIHEDRALS' ) THEN
 120     READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 120 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 120
         END IF
         IF(nword.EQ.0) GOTO 120
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 120
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF

         lptor=lptor+1

c---     Out of boundaries: STOP 
         IF(lptor .GT. nmax_ptor) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='P-Torsion'// err_dim//
     &           ' Change "torsp" in parameters.h'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF
         
c---     Charmm-Like torsion 
         IF(nword .EQ. 6) THEN
            altor(1,lptor)=strngs(1)(1:7)
            altor(2,lptor)=strngs(2)(1:7)
            altor(3,lptor)=strngs(3)(1:7)
            altor(4,lptor)=strngs(4)(1:7)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) ptor1(lptor)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) dummy
            ntor2(lptor)=IDNINT(dummy)
            GOTO 120

c---     Amber-Like torsion 
         ELSE IF(nword .EQ. 7) THEN
            altor(1,lptor)=strngs(1)(1:7)
            altor(2,lptor)=strngs(2)(1:7)
            altor(3,lptor)=strngs(3)(1:7)
            altor(4,lptor)=strngs(4)(1:7)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) v0
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) en
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) gamma
!           transfer of sign to n to handle negative v0. 
            en=DSIGN(en,v0)
            IF(near0(gamma-180.0D0)) THEN
               v0=DABS(v0)
            ELSE IF(near0(gamma)) THEN
               v0=-DABS(v0)
!           if(abs(en).gt.10.and.lt.20) phase il 60 
!           if(abs(en).gt.20) phase is 90 (not used at the moment). 
!           In fptors the formulat cos(x+60)=0.5*cos(x) + sin(60)*sin(x) is
!           used.  Sign of the sin is determined as in fitors.
            ELSE IF (near0(gamma-60.0D0)) THEN
              v0=DABS(v0)
              if(en.gt.0.d0 ) en = en + 10
              if(en.lt.0.d0 ) en = en - 10
            ELSE IF (near0(gamma-90.0D0)) THEN 
              v0=DABS(v0)
              if(en.gt.0.d0) en = en + 20
              if(en.lt.0.d0) en = en - 20
            ELSE  
               j=0
               k=iline
               call int_str (k,fmt,j)
               errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
               call xerror (errmsg,80,1,222) 
               errmsg= 'P-Torsion Phase can be either 0, 60,90 or 180.'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               GO TO 120
            END IF
            ptor1(lptor)=v0
            ntor2(lptor)=IDNINT(en)              
            GOTO 120
            
c---     God-knows-what-like torsion 
         ELSE IF(nword .EQ. 8) THEN
            altor(1,lptor)=strngs(1)(1:7)
            altor(2,lptor)=strngs(2)(1:7)
            altor(3,lptor)=strngs(3)(1:7)
            altor(4,lptor)=strngs(4)(1:7)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) potj(lptor,1)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) potj(lptor,2)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) potj(lptor,3)
            CALL fndfmt(2,strngs(8),fmt)
            READ(strngs(8),fmt,err=20) potj(lptor,4)
            GOTO 120
         ELSE
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '6, 7 or 8'// err_args // 'TORSION'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GO TO 120
         END IF

c=====keyword  TORSION IMPROPER ========================================

      ELSE IF(( strngs(1).EQ. 'TORSION' .AND. strngs(2) .EQ.
     x        'IMPROPER') .OR. strngs(1).EQ. 'IMPROPER'  ) THEN
 130     READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 130 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 130
         END IF
         IF(nword.EQ.0) GOTO 130
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 130
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF
         lpito=lpito+1

c---     Out of boundaries: STOP 
         IF(lpito .GT. nmax_itor) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='I-Torsion'// err_dim//
     &           ' Change "itorp" in parameters.h'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF
         IF(nword .LT. 6.or.nword.GT.8) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '6 or more'// err_args // 'TORSION IMPROPER'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GOTO 130 
c---     Charmm-Like torsion 
         ELSE IF(nword .EQ. 6) THEN
            alito(1,lpito)=strngs(1)(1:7)
            alito(2,lpito)=strngs(2)(1:7)
            alito(3,lpito)=strngs(3)(1:7)
            alito(4,lpito)=strngs(4)(1:7)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pito1(lpito)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) pito2(lpito)
            pito3(lpito)=1.0D0
            GOTO 130

c---     Kept for compatibility with earlier v. ("harmonic" kywd)
         ELSE IF(nword .EQ. 7) THEN
            alito(1,lpito)=strngs(1)(1:7)
            alito(2,lpito)=strngs(2)(1:7)
            alito(3,lpito)=strngs(3)(1:7)
            alito(4,lpito)=strngs(4)(1:7)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pito1(lpito)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) pito2(lpito)
            IF(strngs(7) .EQ. 'harmonic' .OR. strngs(7) .EQ.
     &           'HARMONIC') THEN
               pito3(lpito)=1.0D0
            ELSE IF(strngs(7) .EQ. 'cosine' .OR. strngs(7) .EQ.
     &              'COSINE') THEN
               j=0
               k=iline
               call int_str (k,fmt,j)
               errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
               call xerror (errmsg,80,1,222) 
               errmsg='Invalid keyword "cosine": ' 
     &              //'7 args. only "harmonic" kywd permitted.'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               GO TO 130
            ELSE
               j=0
               k=iline
               call int_str (k,fmt,j)
               errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
               call xerror (errmsg,80,1,222) 
               errmsg='Invalid keyword "'//strngs(7)(1:7)//'"' 
     &              //'7 args. only "harmonic" kywd permitted.'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               GO TO 130
            END IF
            GOTO 130

c----    amber type torsions 
         ELSE IF(nword .EQ. 8) THEN
            IF(strngs(8) .EQ. 'cosine' .OR. strngs(8) .EQ. 'COSINE')
     &           THEN
               alito(1,lpito)=strngs(1)(1:7)
               alito(2,lpito)=strngs(2)(1:7)
               alito(3,lpito)=strngs(3)(1:7)
               alito(4,lpito)=strngs(4)(1:7)
               CALL fndfmt(2,strngs(5),fmt)
               READ(strngs(5),fmt,err=20) v0
               CALL fndfmt(2,strngs(6),fmt)
               READ(strngs(6),fmt,err=20) en
               CALL fndfmt(2,strngs(7),fmt)
               READ(strngs(7),fmt,err=20) gamma
               en=DSIGN(en,v0)
               IF(near0(gamma-180.0D0)) THEN
                  v0=DABS(v0)
               ELSE IF(near0(gamma)) THEN
                  v0=-DABS(v0)
               ELSE
                  j=0
                  k=iline
                  call int_str (k,fmt,j)
                  errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
                  call xerror (errmsg,80,1,222) 
                  errmsg= 'P-Torsion Phase can be either 0 o 180.'
                  call xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
                  GO TO 130 
               END IF
               pito1(lpito)=v0        
               pito2(lpito)=en
               pito3(lpito)=-1.0D0
            ELSE
               j=0
               k=iline
               call int_str (k,fmt,j)
               errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
               call xerror (errmsg,80,1,222) 
               errmsg='Invalid keyword "'//strngs(7)(1:7)//'"' 
     &              //'8 args. only "cosine" kywd permitted.'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               GOTO 130
            END IF
            GOTO 130
         ELSE
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '6, 7 or 8'// err_args // 'IMPROPER TORSION'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GO TO 130
         END IF

c=====keyword  NONBONDED MIXRULE =======================================

      ELSE IF(strngs(1).EQ. 'NONBONDED' .AND. strngs(2) .EQ.
     x        'MIXRULE') THEN
         iz=0
 140     READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 140 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 140
         END IF
         IF(nword.EQ.0) GOTO 140
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 140
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF

         lpnbd=lpnbd+1

c---     Out of boundaries: STOP 
         IF(lpnbd .GT. atop) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='Nbonded types'// err_dim//
     &           ' Change "_TYP_SOLU" in config.H'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF

         IF(nword .LT. 5) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '5, 6 or 7'// err_args // 'NONBONDED MIXRULE'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            GOTO 140
         ELSE IF(nword .EQ. 5) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) mass(lpnbd)
            ejorg(lpnbd)=0.0D0
            sjorg(lpnbd)=0.0D0
            pnbd3(lpnbd)=0.0D0
            pnbd4(lpnbd)=0.0D0
            GOTO 140
         ELSE IF(nword .EQ. 7) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pnbd3(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) sjorg(lpnbd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) ejorg(lpnbd)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) mass(lpnbd)
            GOTO 140
         ELSE IF(nword .EQ. 6) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pnbd3(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pnbd4(lpnbd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) mass(lpnbd)
            GOTO 140
         END IF

c=====keyword  NONBONDED NOMIXRULE======================================

      ELSE IF(strngs(1).EQ. 'NONBONDED' .AND. strngs(2) .EQ.
     x        'NOMIXRULE') THEN
        if(iz.EQ.1.OR.iz.EQ.0) THEN 
          errmsg
     &         ="NOMIXRULES only be used using *one* prm file"
          CALL xerror(errmsg,80,1,30)
          STOP
        end if
c--------read nomixrule parameters  
          iz=1
 150     READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 150 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 150
         END IF
         IF(nword.EQ.0) GOTO 150
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 150
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            DO i=1,lpnbd
               DO j=i,lpnbd
                  ij=j*(j-1)/2+i
c--------------   read "mixed" i-j parameters after END   
 155              READ(kread,'(a78)',END=600) line(1:78)
                  iline = iline + 1
                  IF(jerr.EQ.1) CALL wrenc(kprint,line)
                  IF(line(1:1) .EQ. '#') GOTO 155
                  CALL parse(line,sep,2,comm,strngs,40,nword,
     x                 iret,errmsg)
                  IF(iret.EQ.1) THEN 
                     j1=0
                     k=iline
                     call int_str (k,fmt,j1)
                     errmsg = err_str // fmt(1:j1) // ':'// line(1:48) 
                     call xerror (errmsg,80,1,222) 
                     errmsg='while parsing line: toomany strings'
                     CALL xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                     go to 155
                  END IF
                  IF(nword .LT. 2) THEN
                     k=0
                     l=iline 
                     call int_str (l,fmt,k)
                     errmsg = err_str // fmt(1:k) // ':'// line(1:48) 
                     call xerror (errmsg,80,1,222) 
                     errmsg= '2'// err_args // 'NONBONDED NOMIXRULE'//
     &                    ' after END kywd'
                     call xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                     go to 155 
                  ELSE
                     CALL fndfmt(2,strngs(1),fmt)
                     READ(strngs(1),fmt,err=20) ecc6(ij)
                     CALL fndfmt(2,strngs(2),fmt)
                     READ(strngs(2),fmt,err=20) ecc12(ij)
                  END IF
               END DO
            END DO
            GOTO 100
         END IF
         lpnbd=lpnbd+1
c---     Out of boundaries: STOP 
         IF(lpnbd .GT. atop) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            CALL xerror(errmsg,80,1,222)
            errmsg='Nbonded types'// err_dim//
     &           ' Change "_TYP_SOLU" in config.H'
            CALL xerror(errmsg,80,1,2)
            STOP
         ENDIF

         IF(nword .LT. 5) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '5'// err_args // 
     &           'NONBONDED NOMIXRULE before END kywd'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 150 
         ELSE IF(nword .EQ. 5) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) mass(lpnbd)
            ejorg(lpnbd)=0.0D0
            sjorg(lpnbd)=0.0D0
            pnbd3(lpnbd)=pnbd1(lpnbd)
            pnbd4(lpnbd)=pnbd2(lpnbd)
            GOTO 150
         ELSE IF(nword .EQ. 7) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pnbd3(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) sjorg(lpnbd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) ejorg(lpnbd)
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) mass(lpnbd)
            GOTO 150
         ELSE IF(nword .EQ. 6) THEN
            alnbd(lpnbd)=strngs(1)(1:7)
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) pnbd1(lpnbd)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) pnbd2(lpnbd)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) pnbd3(lpnbd)
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) pnbd4(lpnbd)
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) mass(lpnbd)
            GOTO 150
         END IF

c=====keyword  H-BOND =================================================

      ELSE IF(strngs(1).EQ. 'H-BOND' ) THEN
 160     READ(kread,'(a78)',END=600) line(1:78)
         iline = iline + 1
         IF(jerr.EQ.1) CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 160 
         iret=0
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,
     x        errmsg)
         iret=0
         IF(iret.EQ.1) THEN 
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg='while parsing line: toomany strings'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 160
         END IF
         IF(nword.EQ.0) GOTO 160
         IF(nword.EQ.1 .AND. strngs(1) .EQ. ' ' ) GOTO 160
         IF(nword.EQ.1 .AND. strngs(1) .EQ. 'END' ) THEN
            GOTO 100
         END IF
         lphyd=lphyd+1
         IF(nword .LT. 4) THEN
            j=0
            k=iline
            call int_str (k,fmt,j)
            errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
            call xerror (errmsg,80,1,222) 
            errmsg= '4'// err_args // 'H-BOND'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
            go to 160 
         ELSE
            alhyd(2,lphyd)=strngs(1)(1:7)
            alhyd(1,lphyd)=strngs(2)(1:7)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) phyd1(lphyd)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) phyd2(lphyd)
            GOTO 160
         END IF

c=====Blank Line  =======================================================

      ELSE IF(strngs(1) .EQ. ' ') THEN
         GOTO 100
      ELSE 
         j=0
         k=iline
         call int_str (k,fmt,j)
         errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
         call xerror (errmsg,80,1,222) 
         errmsg= err_unr// strngs(1)(1:20)
         call xerror (errmsg,80,1,2) 
         STOP 
      END IF

 600  CONTINUE

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j)//' ERRORS WHILE EXECUTING RDAPRM'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN
         CALL xerror(errmsg,80,1,30)
         errmsg='MORE THAN 99 ERRORS WHILE EXECUTING RDAPRM'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      iret=0 

      if(ffchar=="R") THEN       ! this means that a PRM file has been already loaded 
! if so prepend "*" to get no overlap between ligand and protein
! PRMs
        do i=loc_pbon+1,lpbon
          albon(1,i)(1:7)="*"//albon(1,i)(1:6)
          albon(2,i)(1:7)="*"//albon(2,i)(1:6)
        end do
        do i=loc_pbnd+1,lpbnd
          albnd(1,i)(1:7)="*"//albnd(1,i)(1:6)
          albnd(2,i)(1:7)="*"//albnd(2,i)(1:6)
          albnd(3,i)(1:7)="*"//albnd(3,i)(1:6)
        end do
        do i=loc_ptor+1,lptor
          if(altor(1,i)(1:1).NE."x")altor(1,i)(1:7)="*"//altor(1,i)(1:6)
          if(altor(2,i)(1:1).NE."x")altor(2,i)(1:7)="*"//altor(2,i)(1:6)
          if(altor(3,i)(1:1).NE."x")altor(3,i)(1:7)="*"//altor(3,i)(1:6)
          if(altor(4,i)(1:1).NE."x")altor(4,i)(1:7)="*"//altor(4,i)(1:6)
        end do
        do i=loc_pito+1,lpito
          if(alito(1,i)(1:1).NE."x")alito(1,i)(1:7)="*"//alito(1,i)(1:6)
          if(alito(2,i)(1:1).NE."x")alito(2,i)(1:7)="*"//alito(2,i)(1:6)
          if(alito(3,i)(1:1).NE."x")alito(3,i)(1:7)="*"//alito(3,i)(1:6)
          if(alito(4,i)(1:1).NE."x")alito(4,i)(1:7)="*"//alito(4,i)(1:6)
        end do
        do i=loc_pnbd+1,lpnbd
          alnbd(i)(1:7)="*"//alnbd(i)(1:6)
        end do
        do i=loc_phyd+1,lphyd
          alhyd(2,i)(1:7)="*"//alhyd(2,i)(1:6)
          alhyd(1,i)(1:7)="*"//alhyd(1,i)(1:6)
        end do
      end if

      RETURN

c==============================================================================
c     Wrong format  - Internal Read error - Tab character
c==============================================================================


 20   CONTINUE
      j=0
      k=iline
      call int_str (k,fmt,j)
      errmsg = err_str // fmt(1:j) // ':'// line(1:48) 
      call xerror (errmsg,80,1,30) 
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)

      STOP
      END
