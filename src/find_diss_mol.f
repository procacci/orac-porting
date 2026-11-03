      SUBROUTINE find_diss_mol(diss_mol,diss_atoms,atoms_diss,mol_diss
     &     ,diss_list,nprot,protl)

************************************************************************
*   Time-stamp: <97/04/28 13:54:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL diss_mol,diss_atoms
      INTEGER atoms_diss(2,2),mol_diss(2),diss_list(2,*),nprot,protl(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,k,m,count

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO n=1,nprot
         m=protl(1+count)
         IF(n .EQ. mol_diss(1)) THEN
            diss_list(1,1)=m
            DO k=1,m
               diss_list(1,k+1)=protl(1+count+k)
            END DO
         ELSE IF(n .EQ. mol_diss(2)) THEN
            diss_list(2,1)=m
            DO k=1,m
               diss_list(2,k+1)=protl(1+count+k)
            END DO
         END IF
         count=count+1+m
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
