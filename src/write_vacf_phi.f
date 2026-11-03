      SUBROUTINE write_vacf_phi(buffer_time,nstep,fstep,vacf,phi)

************************************************************************
*   Time-stamp: <97/12/15 13:52:59 marchi>                             *
*                                                                      *
*  Write velocity autocorrelation function and spectra for solvent     *
*  and solute                                                          *
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

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER buffer_time,nstep
      REAL*8  fstep,vacf(buffer_time,*),phi(buffer_time,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j
      REAL*8  c,nub,tstep,wd
      PARAMETER(c=2.997925D+10)

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(kvaf,100) 
      nub=1.0D0/(DFLOAT(nstep)*2.0D0*fstep*c*1.0D-15)
      DO j=1,2
         IF(DABS(vacf(1,j)) .GT. 1.0D-6) THEN
            IF(j .EQ. 1) WRITE(kvaf,200)
            IF(j .EQ. 2) WRITE(kvaf,300)
            DO i=1,nstep
               tstep=DFLOAT(i-1)*fstep
               wd=DFLOAT(i-1)*nub
               WRITE(kvaf,400) tstep,vacf(i,j),wd,phi(i,j)
            END DO
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(72('*')/'*',70(' '),'*'/
     &     '*         Velocity autocorrelation function and power',
     &     ' spectrum         *'/'*',70(' '),'*'/
     &     '*             Time in femtoseconds and frequency',
     &     ' in cm-1               *'/'*',70(' '),'*'/72('*')/)
200   FORMAT(/'                        > S O L U T E <'/)
300   FORMAT(/'                        > S O L V E N T <'/)
400   FORMAT(3x,f11.3,1x,e15.7,4x,f11.3,1x,e15.7)
      RETURN
      END
