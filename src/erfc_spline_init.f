      SUBROUTINE erfc_spline_init(rspoff,alphal,erfc_bin,mspline
     &     ,erfc_arr,work,iret,errmsg,kcut,rcut_corr,lerf)

************************************************************************
*   Time-stamp: <98/02/24 19:30:23 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Sep  8 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  rspoff,alphal,erfc_bin,erfc_arr(4,*),work(*),kcut
     &     ,rcut_corr
      INTEGER mspline,iret
      CHARACTER*80 errmsg
      logical*4  lerf

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 erftbdns,rcut
      INTEGER mxerftab,imax

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      erftbdns=1.0D0/erfc_bin
      rcut=DSQRT(rspoff)+2.0D0
      if(rcut.lt.rcut_corr) rcut=rcut_corr
      mxerftab=rcut/erfc_bin
      IF(mxerftab .GT. mspline) THEN
         iret=1
         errmsg='In ERFC_SPLINE_INIT: erfc spline dimensions exceed'
     &        //' physical dimensions. Abort.'
         RETURN
      END IF
      CALL fill_erf_table(erftbdns,mxerftab,erfc_arr,work,alphal,kcut
     &     ,rcut,lerf)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
