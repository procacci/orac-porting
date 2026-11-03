      SUBROUTINE read_setup(err_open,err_args,err_end,err_unr
     &           ,err_fnf)

************************************************************************
*   Time-stamp: <98/01/26 16:12:53 marchi>                             *
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
      
      INTEGER iret
      CHARACTER*22 err_open
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,nco,nsevere,nwarning,j
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      CHARACTER*15 err_fnf
      LOGICAL  exist
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nsevere = 0
      nwarning = 0
      line(79:80)='  '

c=======================================================================
c     Environment parser starts here 
c=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)

      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF
      
c==== Command  READ_PDB================================================

      IF(strngs(1).EQ. 'READ_PDB' ) THEN
         pdb_read=.TRUE.
         if (strngs(2).EQ."many") THEN
           errmsg=" ...Expecting multiple PDBs for NE alchemy"
           CALL xerror(errmsg,80,1,11)
         ELSE
           INQUIRE(FILE=strngs(2),EXIST=exist)
           IF(exist) THEN
             pdbfile=strngs(2)
             CALL openf(kconf,strngs(2),'FORMATTED','OLD',0)
           ELSE
             nsevere = nsevere + 1
             errmsg='PDB'// err_fnf
             CALL xerror(errmsg,80,1,30)
           END IF
         ENDIF 

c==== Command  RESET_CM ==============================================

      ELSE IF(strngs(1).EQ. 'RESET_CM' ) THEN
         rescm=.TRUE.
         
c==== Command  SOLUTE==================================================

      ELSE IF(strngs(1) .EQ. 'SOLUTE') THEN
         IF(nword .EQ. 1) THEN
            slt_exist=.TRUE.
         ELSE IF(nword .EQ. 2) THEN
            IF(strngs(2) .EQ. 'ON') THEN
               slt_exist=.TRUE.
            ELSE IF(strngs(2) .EQ. 'OFF') THEN
               slt_exist=.FALSE.
            ELSE 
               nsevere = nsevere + 1
               errmsg= err_args(2)//'2'
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            errmsg=err_unr(3) //strngs(3)(1:8)
            nsevere = nsevere + 1
            CALL xerror(errmsg,80,1,30)
         END IF
               
c==== Command  SOLVENT=================================================

      ELSE IF(strngs(1) .EQ. 'SOLVENT') THEN
         IF(nword .EQ. 1) THEN
            slt_exist=.TRUE.
         ELSE IF(nword .EQ. 2) THEN
            IF(strngs(2) .EQ. 'ON') THEN
               slt_exist=.TRUE.
            ELSE IF(strngs(2) .EQ. 'OFF') THEN
               slt_exist=.FALSE.
            ELSE 
               nsevere = nsevere + 1
               errmsg= err_args(2)//'2'
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            errmsg=err_unr(3) //strngs(3)(1:8)
            nsevere = nsevere + 1
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  REPLICATE===============================================

      ELSE IF(strngs(1).EQ. 'REPLICATE' ) THEN
         replicate=.TRUE.
         IF(nword .EQ. 2) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) icl
            icm=icl
            icn=icl
         ELSE IF(nword .EQ. 4) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) icl
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) icm
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) icn
         ELSE
            nsevere = nsevere + 1
            errmsg= err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  CHANGE_CELL=============================================

      ELSE IF(strngs(1).EQ. 'CHANGE_CELL' ) THEN
         change_cell=.TRUE.

c==== Command  READ_CO=================================================

      ELSE IF(strngs(1).EQ. 'READ_CO' ) THEN
         read_co=.TRUE.
         nco=0
c------- read the line
200      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100

         ELSE IF(nword .EQ. 3) THEN
            nco=nco+1
            if(nco .GT. 3) THEN 
               nsevere = nsevere + 1
               errmsg= 'CO matrix is 3X3. Rows exceed 3. Abort.'
               CALL xerror(errmsg,80,1,30)
            ELSE
               CALL fndfmt(2,strngs(1),fmt)
               READ(strngs(1),fmt,err=20) co(nco,1)
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) co(nco,2)
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) co(nco,3)
            ENDIF    
            GOTO 200
         ELSE
            nsevere = nsevere + 1
            errmsg= 'READ_CO excepts in arguments only three columns. '
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  CRYSTAL=================================================

      ELSE IF(strngs(1).EQ. 'CRYSTAL' ) THEN
         IF(nword .EQ. 2) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) aaxis
            baxis=aaxis
            caxis=aaxis
            alf=90.0D0
            bet=90.0D0
            gam=90.0D0
         ELSE IF(nword .EQ. 4) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) aaxis
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) baxis
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) caxis
            alf=90.0D0
            bet=90.0D0
            gam=90.0D0
         ELSE IF(nword .EQ. 7) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) aaxis
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) baxis
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) caxis
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) alf
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) bet
            CALL fndfmt(2,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) gam
         ELSE IF(nword.eq.1) THEN
            read_remark=.true.
         ELSE
            nsevere = nsevere + 1
            errmsg= err_args(1)// '1' 
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  TEMPLATE================================================

      ELSE IF(strngs(1).EQ. 'TEMPLATE' ) THEN
         if(nword.eq.2) THEN 
            template=.TRUE.
            CALL uscrpl(strngs(2),80)
            INQUIRE(FILE=strngs(2),EXIST=exist)
            IF(exist) THEN
               CALL openf(ktemplate,strngs(2),'FORMATTED','OLD',0)
            ELSE
               nsevere = nsevere + 1
               errmsg= 'Template'//err_fnf
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg= err_args(1)// '1' 
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  RECONSTRUCT=============================================

      ELSE IF(strngs(1).EQ. 'RECONSTRUCT' ) THEN
         recstrc=.TRUE.

c==== Blank Line ======================================================

      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1)// strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== No &END .or. keyword found ======================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

      ELSE
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

      GOTO 100
600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   if syntax errors abort
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j)//' ERRORS WHILE EXECUTING READ_SETUP'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN
         CALL xerror(errmsg,80,1,30)
         errmsg='MORE THAN 99 ERRORS WHILE EXECUTING READ_SETUP'
         call xerror(errmsg,80,1,2)
         STOP
      END IF

      RETURN

c==============================================================================
c     Wrong format  - Internal Read error - Tab character
c==============================================================================


 20   CONTINUE
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
