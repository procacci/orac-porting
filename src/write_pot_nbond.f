      SUBROUTINE write_pot_nbond(ecc6,ecc12,ecc146,ecc1412,nbtype,mass
     &     ,fudge,beta,nato)

************************************************************************
*   Time-stamp: <98/03/07 11:17:14 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Dec 11 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato,ndim
      INTEGER nbtype(*)
      CHARACTER*7 beta(*)
      REAL*8 fudge,ecc6(*),ecc12(*),ecc146(*),ecc1412(*),mass(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER indexa(m1)
      REAL*8  siga(m1),sigb(m1),epsa(m1),epsb(m1)
      COMMON /rag1/ siga,sigb,epsa,epsb,indexa

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ii,natob,ia,ib,j
      REAL*8  fact,c6,c12,aux
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      aux=unite*avogad/1000.0D0/4.184D0
      fact=(2.0D0**(1.0D0/6.0D0))*0.5D0
      DO i=1,nato
         indexa(i)=0
      END DO
      natob=0
      DO i=1,nato
         ii=nbtype(i)*(nbtype(i)+1)/2
         c6=ecc6(ii)
         c12=ecc12(ii)
         IF(c6 .NE. 0.0D0 .AND. c12 .NE. 0) THEN
            siga(i)=(c12/c6)**(1.0D0/6.0D0)*fact
            epsa(i)=0.25D0*(c6*c6/c12)*aux
         ELSE
            siga(i)=0.0D0
            epsa(i)=0.0D0
         END IF

         c6=ecc146(ii)
         c12=ecc1412(ii)
         IF(c6 .NE. 0.0D0 .AND. c12 .NE. 0) THEN
            sigb(i)=(c12/c6)**(1.0D0/6.0D0)*fact
            epsb(i)=0.25D0*(c6*c6/c12)*aux
            IF(DABS(epsb(i)/epsa(i)-fudge) .LT. 1.0D-6 .AND. DABS(sigb(i
     &           )/siga(i)-fudge) .LT. 1.0D-6) THEN
               epsb(i)=0.0D0
               sigb(i)=0.0D0
            END IF
         ELSE
            sigb(i)=0.0D0
            epsb(i)=0.0D0
         END IF
      END DO


      DO i=1,nato
         ia=i
         ok=.TRUE.
         DO j=1,natob
            ib=indexa(j)
            IF(beta(ia) .EQ. beta(ib)) THEN
               ok=.FALSE.
               GOTO 100
            END IF
         END DO
100      CONTINUE
         IF(ok) THEN
            natob=natob+1
            indexa(natob)=ia
         END IF
      END DO

      DO i=1,natob
         ia=indexa(i)
         WRITE(kprint,'(a7,5f11.5)') beta(ia),siga(ia),epsa(ia),sigb(ia)
     &        ,epsb(ia),mass(ia)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
