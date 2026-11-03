      SUBROUTINE open_file_dump(uniit,filename,recla,reclb
     &     ,divide_records,atom_record,nbyte,nrec,occupy,analys,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <97/11/29 18:34:33 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr 21 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER uniit,iret,recla,reclb,divide_records,atom_record,nbyte
      CHARACTER*80 filename,errmsg
      LOGICAL occupy,analys

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER count,i,j,map,recl1,recl2,naux
      INTEGER nword,nsevere,nrec,norec,strblk,nato
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist,ok
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      recl1=recla
      recl2=reclb
      nato=reclb


      WRITE(kprint,10000) 
      
      INQUIRE(FILE=filename,EXIST=exist)
      IF(exist) THEN
         CALL openf(kcnfig,filename,'FORMATTED','OLD',0)
      ELSE
         iret=1
         errmsg=
     &' Input auxiliary file does not exist. Abort.'
         RETURN
      END IF

*=======================================================================
*--- Read the auxiliary file -------------------------------------------      
*=======================================================================

      nsevere=0
      count=0
100   READ(kcnfig,'(a80)',ERR=400,END=200) line
      IF(line(1:1) .EQ. '#') GOTO 100 
      IF(line(1:1) .EQ. ' ') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret .EQ. 1) THEN 
         errmsg='While parsing line: Too many strings'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 200
      END IF
      count=count+1
      file_names(count)=strngs(1)
      IF(nword .EQ. 4) THEN
         CALL fndfmt(1,strngs(2),fmt)
         READ(strngs(2),fmt,err=400) no_records(count)
         CALL fndfmt(1,strngs(3),fmt)
         READ(strngs(3),fmt,err=400) divide_records
         CALL fndfmt(1,strngs(4),fmt)
         READ(strngs(4),fmt,err=400) atom_record
         ok=.TRUE.
      ELSE IF(nword .EQ. 1) THEN
         ok=.FALSE.
      ELSE
         errmsg='While parsing line: Four substrings are expected.'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 200
      END IF
      GOTO 100
200   CONTINUE

      IF(nsevere .GT. 0) THEN
         call int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)/ /' ERRORS WHILE EXECUTING OPEN_FILE_DUMP'
         call xerror(errmsg,80,1,2)
         RETURN
      END IF

      IF(count .EQ. 0) THEN
         iret=1
         errmsg=
     &        'Auxiliary file is empty. Abort'
         RETURN
      ELSE IF(count .EQ. 1) THEN
         iret=1
         errmsg=
     &        'Auxiliary file must contain at least two '/ /
     &        'file names. Abort.'
         RETURN
      END IF

*=======================================================================
*---- Compute record numbers and record length for the trajectory ------
*---- files ------------------------------------------------------------
*=======================================================================
            
      nwrite_dump=count

*---- The file was rewritten by program 

      IF(ok) THEN
         count=0
         DO i=2,nwrite_dump
            count=count+no_records(i)
         END DO
         reclb=atom_record*3*nbyte
      ELSE

*---- The file was written by the user

*---- compute reclb and atom_record and recompute atom_record if 
*---- not a divisor of the number of atoms

         divide_records=reclb/atom_record
         IF(MOD(reclb,atom_record) .NE. 0) THEN
            divide_records=divide_records+1
         END IF
         reclb=nbyte*3*atom_record

         norec=nrec/(count-1)
         norec=norec*divide_records
         no_records(1)=nrec
         DO i=2,nwrite_dump-1
            no_records(i)=norec
         END DO
         no_records(nwrite_dump)=norec+divide_records*MOD(nrec,count-1)
         count=0
         DO i=2,nwrite_dump
            count=count+no_records(i)
         END DO

*---- Rewrite auxiliary files including any user comments

         REWIND kcnfig
         map=0
300      READ(kcnfig,'(a80)',ERR=400,END=200) line
         map=map+1
         IF(line(1:1) .EQ. '#') GOTO 300 
         IF(line(1:1) .EQ. ' ') GOTO 300 

         map=map-1
         REWIND kcnfig
         DO i=1,map
            READ(kcnfig,'(a80)',ERR=400,END=200) line
         END DO

         WRITE(kcnfig,'(''#  Rewritten by Program'')')
         WRITE(kcnfig,*)
         j=strblk(file_names(1),80)
         naux=0
         WRITE(kcnfig,'(a,5x,3i8)') file_names(1)(1:j),no_records(1)
     &        ,naux,naux
         DO i=2,nwrite_dump
            j=strblk(file_names(i),80)
            WRITE(kcnfig,'(a,5x,3i8)') file_names(i)(1:j),no_records(i)
     &           ,divide_records,atom_record
         END DO
      END IF
      CLOSE(kcnfig)

*---  load reclb in recl2

      recl2=reclb

      IF(nsevere .GT. 0) THEN
         call int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)/ /' ERRORS WHILE EXECUTING OPEN_FILE_DUMP'
         RETURN
      END IF

      IF(count .LT. nrec*divide_records) THEN
         iret=1
         errmsg=
     &        'Number of records of the dump files is smaller than'
     &        / /' the length of the simulation.'
         RETURN
      END IF

c$$$      IF(divide_records*atom_record .LT. nato) THEN
c$$$         iret=1
c$$$         errmsg=
c$$$     &        'Number of atoms per record is smaller than'
c$$$     &        / /' the system atoms. Abort'
c$$$         RETURN
c$$$      END IF

*=======================================================================
*--- Open the trajectory files -----------------------------------------
*=======================================================================

      DO i=1,nwrite_dump
         INQUIRE(FILE=file_names(i),EXIST=exist)
         IF(exist .AND. occupy) THEN
            iret=1
            errmsg=
     &' Trajectory files exist and cannot '
     &           / /'be zeroed. *occupy* keyword improperly used.'
            RETURN
         END IF
         IF((.NOT. exist) .AND. analys) THEN
            iret=1
            errmsg=
     &' Cannot analyze trajectories if the trajectory files do not'
     &/ /' exist.'
            RETURN
         END IF

         IF(i .NE. 1) THEN
            IF(exist) THEN
               CALL openf(kwrite_dump(i),file_names(i),'U','OLD',recl2)
            ELSE
               CALL openf(kwrite_dump(i),file_names(i),'U','NEW',recl2)
            END IF
         ELSE 
            IF(exist) THEN
               CALL openf(kwrite_dump(i),file_names(i),'U','OLD',recl1)
            ELSE
               CALL openf(kwrite_dump(i),file_names(i),'U','NEW',recl1)
            END IF
         END IF
      END DO

      RETURN

400   CONTINUE
      iret=1
      errmsg=
     &'Internal reading error: wrong format?? Tab Character'
      CALL xerror(errmsg,80,1,2)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

10000 FORMAT(/'        <------ Opening Trajectory File ------->'/)
      STOP
      END
