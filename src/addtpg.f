      SUBROUTINE addtpg(knlist,kprint,jerr,xbond,xnbond,rbond,
     x           xtor,xntor,rtor,xitor,xnitor,ritor,iret,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine add a new list of bonds, torsions and            *
*     improper torsions. Use this to redefine the topology of          *
*     a residue already defined.                                       *
*                                                                      *
*======================================================================*
*                                                                      *
*     KNLIST   :  Unit from which to read input data            (I)    *
*     KPRINT   :  Unit used to print output                     (I)    *
*     JERR     :  Hide output if JERR .NE. 1                    (I)    *
*     XBOND    :  List of bonds for a residue by atom label     (O)    *
*                 >> CHARACTER*7 XBOND(2,*) <<                         *
*     XNBOND   :  Number of new bonds read in by ADDTPG         (O)    *
*     RBOND    :  List of bonds for a residue by atom number    (O)    *
*                 >> INTEGER RBOND(2,*) <<                             *
*     XTOR     :  List of torsion for a residue by atom label   (O)    *
*                 >> CHARACTER*7 XTOR(4,*) <<                          *
*     XNTOR    :  Number of new torsion read in by ADDTPG       (O)    *
*     RTOR     :  List of bonds for a residue by atom number    (O)    *
*                 >> INTEGER RTOR(4,*) <<                              *
*     XITOR    :  List of torsion for a residue by atom label   (O)    *
*                 >> CHARACTER*7 XITOR(4,*) <<                         *
*     XNITOR   :  Number of new torsion read in by ADDTPG       (O)    *
*     RITOR    :  List of bonds for a residue by atom number    (O)    *
*                 >> INTEGER RITOR(4,*) <<                             *
*                                                                      *
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  25/10/92 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Universite' de Paris-Sud        *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER knlist,kprint,jerr,xnbond,xntor,xnitor,iret,rbond(2,*),
     x        rtor(4,*),ritor(4,*)
      CHARACTER*7  xbond(2,*),xtor(4,*),xitor(4,*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,n,number,number1,number2
      INTEGER nword
      CHARACTER*8 char
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1,fior
      DATA sep/' ',','/comm/'#','!'/
      CHARACTER*80 errmsg

*==================== EXECUTABLE STATEMENTS ============================

      xnbond=0
      xntor=0
      xnitor=0
      line(79:80)='  '
100   READ(knlist,'(a78)',END=600) line(1:78)
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.0) THEN
          IF(strngs(1) .EQ. 'bond') THEN
              xnbond=xnbond+1
              IF(strngs(4) .NE. 'residue') GOTO 1000
              IF(strngs(2)(1:1) .EQ. '1' .AND. strngs(3)(1:1) .EQ.
     x           '1') THEN
                  xbond(1,xnbond)(1:7)=strngs(2)(2:7)
                  xbond(2,xnbond)(1:7)=strngs(3)(2:7)
                  CALL fndfmt(1,strngs(5),fmt)
                  READ(strngs(5),fmt) number
                  rbond(1,xnbond)=number
                  rbond(2,xnbond)=number

              ELSE IF(strngs(2)(1:1) .EQ. '1' .AND.  strngs(3)(1:1) 
     x           .EQ. '2') THEN
                  xbond(1,xnbond)(1:7)=strngs(2)(2:7)
                  xbond(2,xnbond)(1:7)=strngs(3)(2:7)
                  CALL fndfmt(1,strngs(5),fmt)
                  READ(strngs(5),fmt) number
                  rbond(1,xnbond)=number
                  CALL fndfmt(1,strngs(6),fmt)
                  READ(strngs(6),fmt) number
                  rbond(2,xnbond)=number
              ELSE 
                  GOTO 1000
              END IF
              
          ELSE IF(strngs(1) .EQ. 'torsion') THEN
              IF(strngs(6) .NE. 'RESIDUE') GOTO 1000
              xntor=xntor+1
              xtor(1,xntor)(1:7)=strngs(2)(2:7)
              xtor(2,xntor)(1:7)=strngs(3)(2:7)
              xtor(3,xntor)(1:7)=strngs(4)(2:7)
              xtor(4,xntor)(1:7)=strngs(5)(2:7)
              IF(strngs(2)(1:1) .EQ. '1' .AND. 
     x           strngs(3)(1:1) .EQ. '1' .AND. 
     x           strngs(4)(1:1) .EQ. '1' .AND. 
     x           strngs(5)(1:1) .EQ. '1') THEN
                  CALL fndfmt(1,strngs(7),fmt)
                  READ(strngs(7),fmt) number
                  rtor(1,xntor)=number
                  rtor(2,xntor)=number
                  rtor(3,xntor)=number
                  rtor(4,xntor)=number
              ELSE 
                  CALL fndfmt(1,strngs(7),fmt)
                  READ(strngs(7),fmt) number1
                  CALL fndfmt(1,strngs(8),fmt)
                  READ(strngs(8),fmt) number2
                  DO i=1,4
                      IF(strngs(i+1)(1:1) .EQ.'1') THEN
                          rtor(i,xntor)=number1
                      ELSE IF(strngs(i+1)(1:1) .EQ. '2') THEN
                          rtor(i,xntor)=number2
                      ELSE
                          GOTO 1000
                      END IF
                  END DO
              END IF


          ELSE IF(strngs(1) .EQ. 'i_torsion') THEN
              IF(strngs(6) .NE. 'RESIDUE') GOTO 1000
              xnitor=xnitor+1
              xitor(1,xnitor)(1:7)=strngs(2)(2:7)
              xitor(2,xnitor)(1:7)=strngs(3)(2:7)
              xitor(3,xnitor)(1:7)=strngs(4)(2:7)
              xitor(4,xnitor)(1:7)=strngs(5)(2:7)
              IF(strngs(2)(1:1) .EQ. '1' .AND. 
     x           strngs(3)(1:1) .EQ. '1' .AND. 
     x           strngs(4)(1:1) .EQ. '1' .AND. 
     x           strngs(5)(1:1) .EQ. '1') THEN
                  CALL fndfmt(1,strngs(7),fmt)
                  READ(strngs(7),fmt) number
                  ritor(1,xnitor)=number
                  ritor(2,xnitor)=number
                  ritor(3,xnitor)=number
                  ritor(4,xnitor)=number
              ELSE 
                  CALL fndfmt(1,strngs(7),fmt)
                  READ(strngs(7),fmt) number1
                  CALL fndfmt(1,strngs(8),fmt)
                  READ(strngs(8),fmt) number2
                  DO i=1,4
                      IF(strngs(i+1)(1:1) .EQ.'1') THEN
                          ritor(i,xnitor)=number1
                      ELSE IF(strngs(i+1)(1:1) .EQ. '2') THEN
                          ritor(i,xnitor)=number2
                      ELSE
                          GOTO 1000
                      END IF
                  END DO
              END IF
          ELSE IF(strngs(1) .EQ. ' ') THEN
              GOTO 100
              
          ELSE IF(strngs(1) .EQ. 'END') THEN
              RETURN

          ELSE
              GOTO 2000
          ENDIF
          GOTO 100
      END IF

600   iret=1
      errmsg=' Error in ADDTPG: E.O.F. encountered. Abort.'
      RETURN

1000  iret=1
      errmsg=' Error in ADDTPG: directive keywords out of sequence.'//
     x  ' Abort.'

      RETURN

2000  iret=1
      errmsg=' Error in ADDTPG: directive NOT found. Abort.'      

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
