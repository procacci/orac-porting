      SUBROUTINE check_length_sim(start_anl,stop_anl,buffer_time
     &     ,buffer_fft,start_time,end_time,length_run,length_tot
     &     ,length_fft,iret,errmsg)

************************************************************************
*   Time-stamp: <97/12/05 11:27:43 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Dec  3 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER start_anl,stop_anl,buffer_time,length_run,length_tot
     &     ,length_fft,start_time,end_time,iret,buffer_fft
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*8 string

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      start_time=start_anl
      end_time=stop_anl
      length_run=stop_anl
      length_tot=end_time-start_time+1
      IF(MOD(length_tot,2) .NE. 0) THEN
         end_time=end_time-1
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
