      SUBROUTINE read_solute(err_args,err_unr,err_end,err_fnf,err_open)

************************************************************************
*   Time-stamp: <98/03/20 15:21:13 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Nov 18 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*15 err_fnf
      CHARACTER*22 err_open

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,nsevere,j
      INTEGER iret
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      j=0
      nsevere = 0
      line(79:80)='  '
      protei=.TRUE.
      ncofactor = 0
      slt_exist=.TRUE.

c=======================================================================
c     Environment parser starts here 
c=======================================================================

 100  READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Command  SPACE_GROUP=============================================

      IF(strngs(1).EQ. 'SPACE_GROUP' ) THEN
         IF(nword.ge.4) THEN 
            sgroup=.TRUE.
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kgroup,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  nsevere = nsevere + 1
                  errmsg='Space-group-symmetry'//err_fnf
                  CALL xerror(errmsg,80,1,30)
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
            END IF
            CALL get_sgr(line,cgroup)
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1)//'3'
            CALL xerror(errmsg,80,1,30)
         END IF   

c==== Command  EXIST ==================================================

      ELSE IF(strngs(1).EQ. 'EXIST' ) THEN
         
c==== Command  COORDINATES=============================================

      ELSE IF(strngs(1).EQ. 'COORDINATES' ) THEN
         INQUIRE(FILE=strngs(2),EXIST=exist)
         IF(exist) THEN
           pdbfile=strngs(2)
           CALL openf(kcoord_slt,strngs(2),'FORMATTED','OLD',0)
         ELSE
            nsevere = nsevere + 1
            errmsg='PDB'// err_fnf
            CALL xerror(errmsg,80,1,30)
         END IF
         
c==== Command  COMMAND SCALE_CHARGES =================================

      ELSE IF(strngs(1).EQ. 'SCALE_CHARGES' ) THEN
         scharge=.TRUE.
         IF(nword .NE. 1) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nprot_charges
            IF(nword .NE. nprot_charges+2) THEN
               errmsg=err_args(1) // strngs(2) //' after arg. #2' 
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               DO i=1,nprot_charges
                  CALL fndfmt(1,strngs(2+i),fmt)
                  READ(strngs(2+i),fmt,err=20) prot_charges(i)
               END DO
            END IF
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  DEF_SOLUTE=============================================

      ELSE IF(strngs(1).EQ. 'DEF_SOLUTE' ) THEN
         anprot=.TRUE.
         annpro=annpro+1
         IF(nword .NE. 3) THEN
            errmsg=err_args(1)//'2'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            DO i=2,3
               CALL fndfmt(1,strngs(i),fmt)
               READ(strngs(i),fmt,err=20) anpoint(i-1,annpro)
            END DO       
         END IF

c==== Command  DEF_COFACTOR===========================================

      ELSE IF(strngs(1).EQ. 'DEF_COFACTOR' ) THEN
         ncofactor = ncofactor+1
         IF(ncofactor.le.20) THEN
            IF(nword .NE. 3) THEN
               errmsg=err_args(1)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) nsegs(1,ncofactor) 
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) nsegs(2,ncofactor) 
            END IF   
         ELSE
            errmsg='Toomany cofactors defined; Max is 20'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Blank line =====================================================

      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

c==== &END or Bad directive ==========================================

      ELSE
c---     could not fine COMMAND or &END
         errmsg= err_unr(1)// strngs(1)(1:8)// err_end
         call xerror(errmsg,80,1,30)
         nsevere = nsevere + 1 
      END IF

      GOTO 100
 600  CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         call int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE IN COMMAND SOLUTE'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE IN COMMAND SOLUTE'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN

 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? Tab Character??' 
      CALL xerror(errmsg,80,1,2)
      STOP
      END
