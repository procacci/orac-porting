      SUBROUTINE copy_protl(protl,protlb,nprot,nprotb)

************************************************************************
*   Time-stamp: <98/01/21 16:12:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jan 21 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER protl(*),protlb(*),nprot,nprotb

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER count,i,m

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nprotb=nprot
      count=0
      DO i=1,nprot
         m=protl(1+count)
         count=count+1+m
      END DO
      DO i=1,count
         protlb(i)=protl(i)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
