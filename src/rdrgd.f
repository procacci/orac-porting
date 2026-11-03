      SUBROUTINE rdrgd(kprint,knlist,jerr,na,stpr,stsc,nstpr,nstsc,iret,
     x                 iline,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine reads the additional rigid constraints           *
*     associated with a given residue. Can be called many              *
*     times.                                                           *
*                                                                      *
*     KPRINT  :  Unit to which write data                              *
*     KNLIST  :  Unit from which to read data                          *
*     JERR    :  Echo flag. If jerr=1 echo input. Do nothing           *
*                otherwise                                             *
*     NA      :  Number of primary sites on the list                   *
*     STPR    :  List of primary sites atom labels.                    *
*                >> character*7 STPR(*) <<                             *
*     STSC    :  List of primary sites atom labels.                    *
*                >> character*7 STSC(*) <<                             *
*     NSTPR   :  Number of primary sites                               *
*     NSTSC   :  Number of secondary sites                             *
*                                                                      *
*----- Last update  11/10/92 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Uniersite' Peris-Sud            *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER na,kprint,knlist,jerr,iline
      INTEGER nstpr,nstsc
      CHARACTER*7 stpr(*),stsc(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m,nword,iret
      CHARACTER*80 line,strngs(40),errmsg
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1,fior
      DATA sep/' ',','/comm/'#','!'/

*==================== EXECUTABLE STATEMENTS ============================


      m=0
      iret=0
      line(79:80)='  '

      
100   READ(knlist,'(a78)',END=600) line(1:78)
      iline = iline + 1
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret .EQ. 1) RETURN
      IF(strngs(1) .EQ. 'end') THEN
          nstsc=m
          nstpr=na
          RETURN
      END IF
          
      IF(m .EQ. 0) THEN
          DO n=1,na
              stpr(n)=strngs(n)(1:7)
          END DO
          DO n=na+1,nword
              stsc(n-na)=strngs(n)(1:7)
          END DO
          m=m+(nword-na)
      ELSE
          DO n=1,nword
              stsc(n+m)=strngs(n)(1:7)
          END DO
          m=m+nword
      END IF
      GOTO 100

600   CONTINUE
      iret=1
      errmsg='Error in RDRGD: EOR found in input stream, abort'
      

*================= END OF EXECUTABLE STATEMENTS ========================
 
      RETURN
      END
