      SUBROUTINE set_ss_array(ss_point,m1,ss_index,nmol,nato_slv,ntap)

************************************************************************
*   Time-stamp: <97/06/29 13:26:10 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 12 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER m1
      INTEGER nato_slv,ntap,nmol,ss_point(m1,*),ss_index(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nts,n,nato

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nts=nato_slv*nmol
      nato=ntap-nts

      IF(nato .EQ. 0) THEN
         ss_point(1,1)=0
      ELSE
         ss_point(1,1)=nato
         DO n=1,nato
            ss_point(n+1,1)=n
         END DO
      END IF
         
      IF(nts .EQ. 0) THEN
         ss_point(1,2)=0
      ELSE
         ss_point(1,2)=nts
         DO n=1,nts
            ss_point(n+1,2)=n+nato
         END DO
      END IF

      IF(ntap .NE. 0) THEN
         DO n=1,nato
            ss_index(n)=1
         END DO
      END IF
      DO n=1,nts
         ss_index(nato+n)=2
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
