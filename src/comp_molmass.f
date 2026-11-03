      SUBROUTINE comp_molmass(nprot,protl,mass,tmass)

************************************************************************
*   Time-stamp: <97/06/26 12:42:14 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar  7 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nprot,protl(*)
      REAL*8  mass(*),tmass(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,count,m
      REAL*8  mtot

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO j=1,nprot
         mtot=0.0D0
         m=protl(1+count)
         DO i=1,m
            mtot=mtot+mass(protl(1+count+i))
         END DO
         tmass(j)=mtot
         count=count+m+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
