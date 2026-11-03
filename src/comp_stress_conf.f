      SUBROUTINE comp_stress_conf(stressd,stressr,prt,oc,volume,unitp
     &     ,press)

************************************************************************
*   Time-stamp: <98/02/13 20:58:09 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Dec  8 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  stressd(3,3),stressr(3,3),prt(3,3),oc(3,3),volume,unitp
     &     ,press

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k
      REAL*8  pressr,pressd

*----------------------- EXECUTABLE STATEMENTS ------------------------*


!=======================================================================
!----- Calculate the pressure from the stress tensor -------------------
!=======================================================================
         
      press = 0.0D0 
      pressd =0.0D0
      pressr =0.0D0
      DO i=1,3
         press=press+stressd(i,i)+stressr(i,i)
         pressd=pressd+stressd(i,i)
         pressr=pressr+stressr(i,i)
      END DO
      press = press*unitp/(3.0D0*volume)
      pressr = pressr*unitp/(3.0D0*volume)
      pressd = pressd*unitp/(3.0D0*volume)
      
!=======================================================================
!----- Calculate the Parrinello-Rahman tensor --------------------------
!=======================================================================

      DO i=1,3
         DO j=1,3
            prt(i,j)=stressd(i,j)+stressr(i,j)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
