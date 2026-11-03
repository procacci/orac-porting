      SUBROUTINE read_solvent(err_open,err_fnf,err_args,err_unr,err_end)

************************************************************************
*   Time-stamp: <97/06/27 18:06:48 marchi>                             *
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
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*15 err_fnf
      CHARACTER*22 err_open

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,j,nsevere,nwarning
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL exist
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nsevere = 0
      nwarning = 0
      line(79:80)='  '
      slv_exist=.TRUE.
      nform=0
      
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

c==== Command  GENERATE  ===============================================
 
      IF(strngs(1).EQ. 'GENERATE' ) THEN
         slv_generate=.TRUE.
         IF(strngs(2) .EQ. 'RANDOMIZE') THEN
            slv_randomize=.TRUE.
            IF(nword .EQ. 5) THEN
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) icl_slv
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) icm_slv
               CALL fndfmt(1,strngs(5),fmt)
               READ(strngs(5),fmt,err=20) icn_slv
            ELSE
               nsevere = nsevere + 1
               errmsg= err_args(1)//'5'
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE 
            slv_randomize=.FALSE.
            IF(nword .EQ. 4) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) icl_slv
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) icm_slv
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) icn_slv
            ELSE
               nsevere = nsevere + 1
               errmsg= err_args(1)//'4'
               CALL xerror(errmsg,80,1,30)
            END IF
         END IF

c==== Command  INSERT==================================================

      ELSE IF(strngs(1).EQ. 'INSERT' ) THEN
         if(nword.eq.2) THEN 
            linser=.TRUE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) radius
          ELSE
            nsevere = nsevere + 1
            errmsg= err_args(1)// '2'
            CALL xerror(errmsg,80,1,30)
          END IF   


c==== Command  COORDINATES=============================================

      ELSE IF(strngs(1).EQ. 'COORDINATES' ) THEN
         INQUIRE(FILE=strngs(2),EXIST=exist)
         IF(exist) THEN
            CALL openf(kcoord_slv,strngs(2),'FORMATTED','OLD',0)
         ELSE
            nsevere = nsevere + 1
            errmsg='PDB'// err_fnf
            CALL xerror(errmsg,80,1,30)
         END IF
         
c==== Command  CELL====================================================

      ELSE IF(strngs(1).EQ. 'CELL' ) THEN
         IF(strngs(2).EQ. 'SC' ) THEN
            nform=1
            
         ELSE IF(strngs(2).EQ. 'FCC' ) THEN
            nform=4
            rmol(1,2)=0.5d0
            rmol(2,2)=-0.5d0
            rmol(2,3)=0.5d0
            rmol(3,3)=-0.5d0
            rmol(1,4)=-0.5d0
            rmol(3,4)=0.5d0
            
         ELSE IF(strngs(2).EQ. 'BCC' ) THEN
            nform=2
            DO i=1,3
               rmol(i,2)=0.5d0
            END DO
         ELSE
            errmsg=' Lattice type not found. ABORT! '
            nsevere = nsevere+ 1
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command  READ_SOLVENT ===========================================

      ELSE IF(strngs(1) .EQ. 'READ_SOLVENT') THEN
         slv_read=.TRUE.
         IF(nword .LT. 2) THEN
            errmsg=err_args(1)//'1'
            nsevere = nsevere + 1
            CALL xerror(errmsg,80,1,30)
         ELSE IF(nword .EQ. 2) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nmol_create
         END IF

c==== Command  REDEFINE ===============================================

      ELSE IF(strngs(1) .EQ. 'REDEFINE') THEN
         slv_redef=.TRUE.
         IF(nword .LT. 3) THEN
            errmsg=err_args(1)//'2'
            nsevere = nsevere + 1
            CALL xerror(errmsg,80,1,30)
         ELSE IF(nword .EQ. 3) THEN
            IF(strngs(2) .EQ. 'TYPEWISE') THEN
               slv_redef_type=strngs(3)(1:8)
c           ELSE IF(strgns(2) .EQ. 'MOLECULEWISE') THEN
            ELSE
               errmsg= err_unr(3) // strngs(2)(1:12)
               nsevere = nsevere + 1
               CALL xerror(errmsg,80,1,30)
            END IF
         END IF


c==== Command  ADD_UNITS =============================================

      ELSE IF(strngs(1) .EQ. 'ADD_UNITS') THEN
         slv_add=.TRUE.
         IF(nword .LT. 2) THEN
            errmsg=err_args(1)//'1'
            nsevere = nsevere + 1
            CALL xerror(errmsg,80,1,30)
         ELSE IF(nword .EQ. 2) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nmol_create
         END IF

c==== Blank Line ======================================================

      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
c==== No &END .or. keyword found ======================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
      ELSE
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF
      GOTO 100

600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE IN COMMAND SOLVENT'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN
         CALL xerror(errmsg,80,1,30)
         errmsg='MORE THAN 99 ERRORS WHILE IN COMMAND SOLVENT'
         call xerror(errmsg,80,1,2)
         STOP
      END IF

c==============================================================================
c     input verification part start here 
c==============================================================================

C---  Check if number of molecules have been provided -------------------------

      IF(slv_read .AND. nmol_create .EQ. 0) THEN
         errmsg='Must specify number of molecules when'//
     &        ' reading coordinates from PDB file'
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

C---  Check if number of molecules have been provided -------------------------

      IF(slv_add .AND. nmol_create .EQ. 0) THEN
         errmsg='Must specify number of molecules when'//
     &        ' reading coordinates from PDB file'
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

C---  Check if CELL has been called at the same time as GENERATE --------------

      IF((slv_generate .AND. nform .EQ. 0) .OR. ( (.NOT. slv_generate)
     &     .AND. nform.NE. 0)) THEN
         errmsg='>GENERATE< and >CELL< need each other:'
     &//' Define the cell type and the no. of replicas'
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

C---  Check if number of molecules have been provided -------------------------

      IF(slv_read .OR. slv_generate) slv_create=.TRUE.
      IF(slv_read .AND. slv_generate) THEN
         errmsg='Cannot read and generate solvent coordinates'
     &//' in the same run.'
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

C---  Check if coordinates have been given ----------------------------------

      IF(kcoord_slv .EQ. 0 .AND. slv_create) THEN
         errmsg='Use command >COORDINATES< to provide template'
     &//' coordinates of the solvent.'
         nsevere = nsevere + 1
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(.NOT. slv_read .AND. .NOT. slv_add) nmol_create=nform*icl_slv
     &     *icm_slv*icn_slv

c--   abort AFTER verifying input 

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE IN COMMAND SOLVENT'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN
         errmsg='MORE THAN 99 ERRORS WHILE IN COMMAND SOLVENT'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
 
      RETURN

c==============================================================================
c     Errors were found
c==============================================================================

 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format???' 
      CALL xerror(errmsg,80,1,2)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
