      SUBROUTINE bend(kprint,knlist,jerr,alpha,n1,nalpha,iret,iline
     &     ,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     All the bends for a given residue are read and a list            *
*     is created. The physical dimension of the list is                *
*     the only argument provided in input.                             *
*                                                                      *
*     ALPHA   :  List of bends for the residue.                        *
*                >> character*7 ALPHA(2,n1) <<                         *
*     N1      :  Physical column dimension of ALPHA.                   *
*     NALPHA  :  Number of bends for the residue.                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update 09/04/89 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL  Free format input package.                             *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER kprint,knlist,jerr,n1,nalpha,iret,iline
      CHARACTER*7 alpha(3,n1)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,i,m
      CHARACTER*8 label,space
      INTEGER nword
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1
      DATA sep/' ',','/comm/'#','!'/

*==================== EXECUTABLE STATEMENTS ============================

      n=0
      line(79:80)='  '

100   READ(knlist,'(a78)') line(1:78)
      iline = iline + 1
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,30,nword,iret,errmsg)
      IF(strngs(1).NE. 'end' ) THEN
	  IF(strngs(1).EQ. ' ') GOTO 100
	  m=nword/3
	  IF(m*3.EQ.nword) THEN
	      DO 10 i=1,nword,3
		  n=n+1
                  alpha(1,n)=strngs(i)(1:7)
                  alpha(2,n)=strngs(i+1)(1:7)
                  alpha(3,n)=strngs(i+2)(1:7)
10            CONTINUE
          ELSE
	      iret=1
	      errmsg=' In BENDS: Needs three atoms to define a bend '
	      RETURN
          END IF
	  GOTO 100
      END IF
      nalpha=n

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
