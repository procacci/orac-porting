      SUBROUTINE comp_nmol(ss_index,nprot,protl,nmol_slt,nmol_slv
     &     ,nato_slt,nato_slv)

************************************************************************
*   Time-stamp: <98/03/17 18:15:22 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Mar 17 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprot,nmol_slt,nmol_slv,nato_slt,nato_slv,protl(*)
     &     ,ss_index(*)
      
*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER count,n,m,itype,sum(2),sum2(2)

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO n=1,2
         sum(n)=0
         sum2(n)=0
      END DO
      DO n=1,nprot
         m=protl(count+1)
         itype=ss_index(protl(count+1+1))
         sum(itype)=sum(itype)+1
         sum2(itype)=sum2(itype)+m
         count=count+m+1
      END DO

      nmol_slt=sum(1)
      nmol_slv=sum(2)
      nato_slt=sum2(1)
      nato_slv=sum2(2)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
