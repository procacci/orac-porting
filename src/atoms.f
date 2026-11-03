      SUBROUTINE atoms(kprint,knlist,jerr,alpha,type,charg1,jgrppt,
     x                 jngrp,n1,nato,iline,iret,errmsg,mass1)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine reads the atom label, types and charges.         *
*     of a given residue. The physical dimension of the lists          *
*     created in output is the only input parameter.                   *
*                                                                      *
*     ALPHA   :  List of atom labels.                                  *
*                >> character*7 ALPHA(N1) <<                           *
*     TYPE    :  List of atom types.                                   *
*                >> character*7 TYPE(N1) <<                            *
*     CHARG1  :  List of atomic charges (used for intramolecular       *
*                interactions).                                        *
*                >> real*8 CHARG1(N1) <<                               *
*     NTPE    :  List of atom types for the intermolecular             *
*                interactions.                                         *
*                >> integer NTPE(N1) <<                                *
*     CHARG2  :  List of atomic charges (used for intermolecular       *
*                interactions).                                        *
*                >> real*8 CHARG2(N1) <<                               *
*     N1      :  Physical dimension of the above arrays.               *
*     NATO    :  Number of atoms for the residue.                      *
*                                                                      *
*----- Last update  09/03/89 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS Free format input package.                             *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n1,nato,kprint,knlist,jerr,jngrp,iline,iret
      INTEGER jgrppt(2,*)
      CHARACTER*7 alpha(n1),type(n1)
      CHARACTER*80 errmsg
      REAL*8 charg1(n1),mass1(n1)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m,nword,j
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1,fior
      DATA sep/' ',','/comm/'#','!'/

*==================== EXECUTABLE STATEMENTS ============================


      n=0
      m=0
      line(79:80)='  '
100   READ(knlist,'(a78)',END=600) line(1:78)
      iline = iline + 1
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg= 'Toomany strings.'
         RETURN
      END IF

      IF(n .EQ. 0 .AND. m .EQ. 0) THEN
          IF(strngs(1) .NE. 'group') THEN
              errmsg='Expected "group" directive' 
              iret = 1
              RETURN
          END IF
      END IF
      IF(strngs(1) .EQ. 'group') THEN
          m=m+1
          IF(n .NE. 0) THEN
              jgrppt(1,m)=n+1
              jgrppt(2,m-1)=n
          ELSE 
              jgrppt(1,m)=n+1
          END IF
          GOTO 100

      ELSE IF(strngs(1).EQ. ' ' ) THEN
          GOTO 100

      ELSE IF(strngs(1) .EQ. 'end') THEN
          jgrppt(2,m)=n

      ELSE 
          n=n+1
          IF(n .GT. n1) THEN
             iret  = 1
             errmsg = 'Toomany atoms defined; dimensions exceed' 
             RETURN
          END IF 
          alpha(n)(1:7)=strngs(1)(1:7)
          type(n)(1:7)=strngs(2)(1:7)
          CALL fndfmt(2,strngs(3),fmt)
          READ(strngs(3),fmt) charg1(n)
          if(nword.eq.4)  THEN 
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt) mass1(n)
          ELSE
            mass1(n) =0.d0
          END IF
	  GOTO 100
      END IF
      nato=n
      jngrp=m

      RETURN

 600  errmsg = 'Unexpected End Of File'
      iret=1 

      RETURN 
*================= END OF EXECUTABLE STATEMENTS ========================

      END
