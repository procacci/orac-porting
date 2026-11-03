      SUBROUTINE write_diffusion(buffer_time,nstep,fstep,diff,nato1
     &     ,nato2)

************************************************************************
*   Time-stamp: <97/12/05 18:55:35 marchi>                             *
*                                                                      *
*  Write root mean square displacement for solvent and solute          *
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

      INTEGER buffer_time,nstep,nato1,nato2
      REAL*8  fstep,diff(buffer_time,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j
      REAL*8  aux,tstep,weight(2)

*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      weight(1)=DFLOAT(nato1)
      weight(2)=DFLOAT(nato2)

      WRITE(kdiff,100)
      DO j=1,2
         IF(weight(j) .GT. 0.0D0) THEN
            aux=1.0D0/weight(j)

            IF(diff(1,j) .LT. 0) diff(1,j)=-diff(1,j)
            IF(j .EQ. 1) WRITE(kdiff,200)
            IF(j .EQ. 2) WRITE(kdiff,300)
            DO i=1,nstep-1
               tstep=DFLOAT(i-1)*fstep
               WRITE(kdiff,400) tstep,DSQRT(diff(i,j)*aux)
            END DO
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(72('*')/'*',70(' '),'*'/
     &     '*                      Root mean square displacement ',
     &     '                  *'/'*',70(' '),'*'/
     &     '*                        Time is in femtoseconds',
     &     '                       *'/'*',70(' '),'*'/72('*')/)
200   FORMAT(/'                        > S O L U T E <'/)
300   FORMAT(/'                        > S O L V E N T <'/)
400   FORMAT(3x,f11.3,1x,e15.7)
      RETURN
      END
