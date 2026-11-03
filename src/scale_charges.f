      SUBROUTINE scale_charges(kprint,nprot_charges,prot_charges
     &     ,chrge,protl,nprot,scharge)

************************************************************************
*   Time-stamp: <97/02/12 10:42:43 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jan 30 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint,nprot,nprot_charges,protl(*),prot_charges(*)
      REAL*8  chrge(*)
      LOGICAL scharge

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,j1,j2,ncount,noff,natomm
      REAL*8  sum,sum1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      natomm=0
      ncount=0
      sum1=0.0D0
      DO i=1,nprot
         noff=protl(ncount+1)
         sum=0.0D0
         DO j=ncount+2,ncount+1+noff
            j1=protl(j)
            sum=sum+chrge(j1)
         END DO
         IF(scharge) THEN
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  sum1=sum1+sum
                  natomm=natomm+noff
               END IF
            END DO
         END IF
         IF(DABS(sum) .GT. 1.0D-6) THEN
            write(kprint,'(5x,''Found Charge'',f10.5,
     &           '' on Solute Molecule '',i6)') sum,i
         END IF
         ncount=ncount+noff+1
      END DO
      IF(scharge) THEN
         sum1=sum1/DFLOAT(natomm)
         ncount=0
         DO i=1,nprot
            noff=protl(ncount+1)
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  DO j1=ncount+2,ncount+1+noff
                     j2=protl(j1)
                     chrge(j2)=chrge(j2)-sum1
                  END DO
               END IF
            END DO
            ncount=ncount+noff+1
         END DO

        
         write(kprint,'(5x,a)')
     &        'Charge Scaling Done               ---->'
         ncount=0
         sum1=0.0D0
         DO i=1,nprot
            noff=protl(ncount+1)
            sum=0.0D0
            DO j=ncount+2,ncount+1+noff
               j1=protl(j)
               sum=sum+chrge(j1)
            END DO
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  sum1=sum1+sum
                  natomm=natomm+noff
               END IF
            END DO
            IF(DABS(sum) .GT. 1.0D-6) THEN
              write(kprint,'(5x,''Found New Charge'',f10.5,
     &              '' on Solute Molecule '',i6)') sum,i
            END IF
            ncount=ncount+noff+1
         END DO

      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
