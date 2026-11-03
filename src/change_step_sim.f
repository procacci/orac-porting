      SUBROUTINE change_step_sim(buffer_time,buffer_fft,fstep,length_tot
     &     ,length_fft,divide,iret,errmsg)

************************************************************************
*   Time-stamp: <97/12/04 13:24:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Dec  4 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER length_tot,length_fft,divide,iret,buffer_time,buffer_fft
      REAL*8  fstep
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*8  string

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      fstep=fstep/DFLOAT(divide)
      length_tot=length_tot*divide
      IF(MOD(length_tot,2) .NE. 0) THEN
         length_tot=length_tot-1
      END IF
      length_fft=length_tot*2
      WRITE(string,'(i8)') length_tot
      IF(length_tot .GT. buffer_time) THEN
         iret=1
         errmsg=
     &        'Allocated space exceeded.'
     &     / /' Increase _BUFFER_TIME_ to at least:'
     &     / /string
      END IF
      WRITE(string,'(i8)') length_fft
      IF(length_fft .GT. buffer_fft) THEN
         iret=1
         errmsg=
     &        'Allocated space exceeded.'
     &     / /' Increase _BUFFER_FFT_ to at least:'
     &     / /string
      END IF
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
