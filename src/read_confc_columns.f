      SUBROUTINE read_confc_columns(xp0,yp0,zp0,xpb,ypb,zpb,nato,iatom
     &     ,length_time,start_time,end_time,divide_records,atom_record
     &     ,iret,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     CONF will dump the hystory of the MD run. The writing on         *
*     disk is performed by the IOLIB routines.                         *
*                                                                      *
*     XP0     :  Solute molecules site coordinates, packed        (I)  *
*     YP0        by molecule.                                          *
*     ZP0        >> real*8  XP0(*), ... <<                             *
*                                                                      *
*     NSTEP   :  Number of the current configuration being        (I)  *
*                dumped.                                               *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS IOLIB                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,divide_records,atom_record,length_time,start_time
     &     ,end_time,iret,iatom
      REAL*8  xp0(*),yp0(*),zp0(*)
      REAL*4  xpb(length_time,*),ypb(length_time,*),zpb(length_time,*)
      CHARACTER*80 errmsg
      
*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,rec,count,uniit,reca,nstart,nend,bookmark,naux,n,m
     &     ,iatomb
      LOGICAL read_file
      SAVE bookmark
      DATA bookmark/-1/

*==================== EXECUTABLE STATEMENTS ============================

      IF(iatom .EQ. 0 .OR. iatom .GT. nato) THEN
         iret=1
         errmsg=
     &        'While reading by column, the atom number required '/ /
     &        'is zero or larger than total.'
      END IF
      IF(start_time .EQ. 0) THEN
         iret=1
         errmsg=
     &        'While reading by column, demanded beginning time '/ /
     &        'step cannot be zero. Abort.'
      END IF
      IF(start_time .GT. length_time) THEN
         iret=1
         errmsg=
     &        'While reading by column, demanded beginning time '/ /
     &        'step is larger than trajectory. Abort.'
      END IF
      IF(end_time .EQ. 0) THEN
         iret=1
         errmsg=
     &        'While reading by column, demanded end time '/ /
     &        'step cannot be zero. Abort.'
      END IF
      IF(end_time .GT. length_time) THEN
         iret=1
         errmsg=
     &        'While reading by column, demanded end time '/ /
     &        'step is larger than trajectory. Abort.'
      END IF

*=======================================================================
*----- Determines which record should be read at any given time --------
*=======================================================================

      naux=MOD(iatom,atom_record)
      IF(naux .EQ. 0) THEN
         naux=iatom/atom_record
      ELSE
         naux=iatom/atom_record+1
      END IF

*----- Determines the atom index referred to the beginning of record ---

      iatomb=iatom-(naux-1)*atom_record

*----- See if it needs to read or data are already in buffer -----------

      read_file=.NOT. (iatom .GE. bookmark .AND. iatom .LE. bookmark
     &     +atom_record) .OR. bookmark .EQ. -1

*=======================================================================
*----- Read file if necessary ------------------------------------------
*=======================================================================

      IF(read_file) THEN
         DO n=start_time,end_time
            rec=divide_records*(n-1)
            rec=rec+naux
            count=0
            i=1
            DO WHILE(count .LT. rec)
               i=i+1
               count=count+no_records(i)
            END DO
            
            uniit=kwrite_dump(i)
            
            DO j=2,i-1
               rec=rec-no_records(j)
            END DO
            READ(unit=uniit,rec=rec,ERR=100) (xpb(n,j),ypb(n,j),zpb(n
     &           ,j),j=1,atom_record)
         END DO
         bookmark=(naux-1)*atom_record
      END IF

*=======================================================================
*----- Copy coordinates to output arrays -------------------------------
*=======================================================================

      DO n=start_time,end_time
         xp0(n)=xpb(n,iatomb)
         yp0(n)=ypb(n,iatomb)
         zp0(n)=zpb(n,iatomb)
      END DO
      RETURN

*================= END OF EXECUTABLE STATEMENTS ========================

100   CONTINUE 
      iret=1
      errmsg=
     &        'While reading by columns, reading error occurred. '/ /
     &        'End of file??'

      RETURN
      END
