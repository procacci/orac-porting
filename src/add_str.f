      SUBROUTINE add_str(in,maxd,add,n,out)

************************************************************************
*   Time-stamp: <95/08/19 13:06:23 marchi>                             *
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

      CHARACTER*1 in(*),add(*),out(*)
      INTEGER n,maxd

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ia,i,strblk     

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ia=strblk(in,maxd)
      DO i=1,ia-1
         out(i)=in(i)
      END DO
      DO i=ia,ia+n-1
         out(i)=add(i-ia+1)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
