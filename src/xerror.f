      SUBROUTINE xerror(errmsg,length,idummy,err)

************************************************************************
*   Time-stamp: <97/12/03 14:29:58 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Mar 27 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER length,err,idummy
      CHARACTER*1  errmsg(length)

#ifdef _MPI_
      include 'mpif.h'
      INTEGER ierr
#endif
*------------------------- LOCAL VARIABLES ----------------------------*
      
      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(err .EQ. 2) THEN
         WRITE(kprint,2) (errmsg(i),i=1,length)
      ELSE IF (err.eq.20) THEN 
         WRITE(kprint,20) (errmsg(i),i=1,length)

      ELSE IF (err.eq.21) THEN 
         WRITE(kprint,21) (errmsg(i),i=1,length)

      ELSE IF (err.eq.25) THEN 
         WRITE(kprint,25) (errmsg(i),i=1,length)
      ELSE IF (err.eq.26) THEN 
         WRITE(kprint,26) (errmsg(i),i=1,length)

      ELSE IF (err.eq.30) THEN 
         WRITE(kprint,30) (errmsg(i),i=1,length)

      ELSE IF (err.eq.11) THEN 
         WRITE(kprint,11) (errmsg(i),i=1,length)

      ELSE IF (err.eq.221) THEN 
         WRITE(kprint,221)

      ELSE IF (err.eq.222) THEN 
         WRITE(kprint,222) (errmsg(i),i=1,length)

      ELSE IF (err.eq.223) THEN 
         WRITE(kprint,223) 

      ELSE    
         WRITE(kprint,1000) (errmsg(i),i=1,length)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

c2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
c=======================================================================
c====  FORMAT WITH ERR = 2: FATAL; to be used just before stopping =====
c====  ALSO RUNTIME ERROR MUST HAVE THIS FORM
c=======================================================================
 

 2    FORMAT(//' ************  T U R T L E   P O W E R !! ********  '/,
     &        ' * ',80A,/ 
     &        ' ************** Fatal Error Program Stops ********   '//)


c=======================================================================


c20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
c=======================================================================
c====  FORMAT WITH ERR = 20: FATAL; program continues =================
c====  TO BE USED AFTER INPUT HAS BEEN READ IN COMPLETELY IN THE VERI- 
c====  FICATION STEP (E.G. VERIFY_INPUT or in the verif part of READ_X 
c=======================================================================

 20   FORMAT ( / '  * * * E R R O R  * * * '/
     &     ,5x,80a)
25    FORMAT ( '** Action: ',80a)
26    FORMAT ( '**        ',80a)
 21   FORMAT ( / '  * * * W A R N I N G * * * '/
     &     ,5x,80a)

c=======================================================================



c30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 20
c=======================================================================
c====  FORMAT WITH ERR = 30: FATAL; program continues =================
c====  TO BE USED WHILE PARSING WHEN CHECKING SYNTAX IN READ_XX
c=======================================================================

 30   FORMAT ('*** F@#!@- ERR:',1x,80a/)
 11   FORMAT ('??? WARNING: ',1x,80a/)
 

c=======================================================================


c1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
c=======================================================================
c===          FORMAT FOR WARNINGS 
c=======================================================================


 1000 FORMAT(
     &     //' ********* W A R N I N G  W A R N I N G ************'/,
     &     ' * ',80A,
     &     /' ******** Recoverable Error Program Continues ******'//)

c=======================================================================


c====222 222 222 222 222 222 222 222 222 222 222 222 222 222 222 222 222
c===          WRITE PLANE STRING errmsg 
c=======================================================================

 221  FORMAT(//' ************  T U R T L E   P O W E R !! ********  ')
 222  FORMAT (' * ',80a,' *')
 223  FORMAT(' ************** Fatal Error Program Stops ********   '//)

101   CONTINUE
      if(err.EQ.2) THEN 
        STOP
#ifdef _MPI_
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr) 
#endif      
      END IF
      RETURN
      END
