      SUBROUTINE get_spectra_vacf(buffer_time,nstep,nmax,vacf,phi,wsave)

************************************************************************
*   Time-stamp: <97/12/04 17:14:53 marchi>                             *
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

      INTEGER buffer_time,nstep,nmax
      REAL*8  vacf(buffer_time,*),phi(buffer_time,*),wsave(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j
      REAL*8  aux

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL dcosti(nstep,wsave)

      DO j=1,2
         aux=1.0D0/vacf(1,j)
         DO i=1,nstep
            vacf(i,j)=vacf(i,j)*aux
         END DO
         DO i=1,nstep
            phi(i,j)=vacf(i,j)
         END DO
         CALL windows(phi(1,j),nstep,nmax)
         CALL dcost(nstep,phi(1,j),wsave)
         DO i=1,nstep
            phi(i,j)=phi(i,j)/DFLOAT(2*nstep-2)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
