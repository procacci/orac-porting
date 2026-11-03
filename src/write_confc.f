      SUBROUTINE write_confc(co,xp0,yp0,zp0,nato,fstep,ninner,nconf
     &     ,divide_records,atom_record)

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

      INTEGER ninner,nconf,nato,divide_records,atom_record
      REAL*8  xp0(*),yp0(*),zp0(*),fstep
      REAL*8  co(3,3)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,rec,count,uniit,nstart,nend,reca

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*----- Write box and timestep ------------------------------------------
*=======================================================================

      rec=ninner/nconf
      uniit=kwrite_dump(1)

      WRITE(unit=uniit,rec=rec) SNGL(fstep),((co(i,j),j=1,3),i=1,3)
      

*=======================================================================
*----- Determines which file to write to -------------------------------
*=======================================================================

      rec=divide_records*(ninner/nconf-1)+1
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

*=======================================================================
*----- Dump the coordinates of the solute molecules when required ------
*=======================================================================

      DO i=1,divide_records
         nstart=(i-1)*atom_record+1
         nend=nstart-1+atom_record
         IF(nend .GT. nato) nend=nato
         reca=rec+i-1
         WRITE(unit=uniit,rec=reca) (SNGL(xp0(j)),SNGL(yp0(j))
     &        ,SNGL(zp0(j)),j=nstart,nend)
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

