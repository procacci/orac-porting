      SUBROUTINE verify_variables(iret,errmsg)

************************************************************************
*   Time-stamp: <96/06/20 19:24:48 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jun 20 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nmax

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nmax=4000000

      IF(ngrp.GT.nmax) THEN 
         WRITE(kprint,100) ngrp
         STOP
      END IF   
      IF(nmol.GT.nmax) THEN 
         WRITE(kprint,200) nmol
         STOP
      END IF   
      
      IF(ngrp.GT.nmax) THEN 
         WRITE(kprint,100) ngrp
         STOP
      END IF   
      IF(nmol.GT.nmax) THEN 
         WRITE(kprint,200) nmol
         STOP
      END IF   
      

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   format(//'in MTSMD: ******** DISASTER ERROR ************'
     &    /'number of group too large. Cannot address integer',I7,'.',
     &    /'action: Change neighbor list in mtsmd to integer*4 or run',
     &    /'        with no RESPA') 
200   format(//'in MTSMD: ******** DISASTER ERROR ************'
     &    /'number of solvent molecules too large. Cannot address',
     &     ' integer',I7,'.',
     &    /'action: Change neighbor list in mtsmd to integer*4 or run',
     &    /'        with no RESPA') 
      RETURN
      END
