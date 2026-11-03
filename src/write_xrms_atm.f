      SUBROUTINE write_xrms_atm(kprint,nato,label,drp,nstep,fstep,ngrp
     &     ,grppt,res)

************************************************************************
*   Time-stamp: <95/12/16 12:08:43 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Aug 19 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER kprint,nato,nstep,ngrp,grppt(2,*),res(*)
      REAL*8  fstep,drp(*)
      CHARACTER*2  label

*------------------------- LOCAL VARIABLES ----------------------------*
      
      REAL*8  fnstep,aux
      INTEGER j,i
      LOGICAL near0
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(kprint,'('' Tstep = '',f12.2,2x,a2)') fstep,label
      fnstep=DFLOAT(nstep)
      DO j=1,ngrp
         DO i=grppt(1,j),grppt(2,j)
            aux=drp(i)/fnstep
            aux=DSQRT(aux)
            IF(.NOT. near0(aux)) THEN
               WRITE(kprint,'(i6,i6,i5,f13.5)') i,j,res(i),aux
            END IF
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
