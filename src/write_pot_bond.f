      SUBROUTINE write_pot_bond(nbond,na,bond,nb,nc,pbond,beta,flag)

************************************************************************
*   Time-stamp: <98/02/24 16:16:46 marchi>                             *
*                                                                      *
*   Write parameter list only for bending and stretching               *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Dec 10 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nbond,na,nb,nc,flag
      INTEGER bond(na,*)
      REAL*8  pbond(nb,*)
      CHARACTER*7 beta(*)
      CHARACTER*80 label

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER work(m3),work1(m3),index(m3)
      REAL*8  pot(m3)
      LOGICAL xx(m3)
      COMMON /rag1/ work,work1,index,pot,xx

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,count,ia,ib,la,lb,lc,ld,nrbond,nl,strblk,count1,ntph
     &     ,tmp,i1,nrrbond
      REAL*8  a,b,c,d,diff,eps,aux,aux2,sig
      CHARACTER*7 l1(4),l2(4)
      LOGICAL ok
      DATA eps/1.0D-6/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(flag .GT. 0) THEN
         aux=unite*avogad/1000.0D0/4.184D0
         aux2=180.0D0/pi
      ELSE
         aux=1.0D0
         aux2=180.0D0/pi
      END IF
      nrbond=0
      DO i=1,nbond
         work(i)=0
      END DO

      DO i=1,nbond
         ok=.TRUE.
         DO j=1,nrbond
            ib=work(j)

*---- Check if parameters are identical 

            count=0
            DO k=1,nc
               a=pbond(i,k)
               b=pbond(ib,k)
               diff=DABS(a-b)
               IF(a .NE. 0.0D0) diff=DABS(a-b)/DABS(a)
               IF(diff .LT. eps) count=count+1
            END DO
            IF(count .EQ. nc) THEN

               DO k=1,na
                  l1(k)=beta(bond(k,i))
                  l2(k)=beta(bond(k,ib))
               END DO 

*---- Check if label sequence is the same

               count1=0
               DO k=1,na
                  IF(l1(k) .EQ. l2(k)) count1=count1+1
               END DO
               
*---- if this is not the case compare it backward
               
               IF(count1 .NE. na) THEN
                  count1=0
                  DO k=1,na
                     IF(l1(na-k+1) .EQ. l2(k)) count1=count1+1
                  END DO
               END IF

*---- If parameters and sequence is identical the parameters are 
*---- already on the list

               IF(count1 .EQ. na) THEN
                  ok=.FALSE.
                  GOTO 100
               END IF
            END IF
         END DO
100      CONTINUE

*---- Accumulate potentials in work

         IF(ok) THEN
            nrbond=nrbond+1
            work(nrbond)=i
         END IF
      END DO
*---- Only for torsions: Decide which potentials have identical type
*---- for 1 and 4 atoms

      IF(nc .EQ. 2 .AND. na .EQ. 4) THEN
         DO i=1,nrbond
            work1(i)=0
            xx(i)=.TRUE.
         END DO
         nrrbond=0
         DO i=1,nrbond
            ia=work(i)
            ok=.TRUE.
            DO j=1,nrrbond
               ib=work1(j)

*---- Check if parameters are identical 

               count=0
               DO k=1,nc
                  a=pbond(ia,k)*aux
                  b=pbond(ib,k)*aux
                  diff=DABS(a-b)
                  IF(a .NE. 0.0D0)  diff=DABS(a-b)/DABS(a)
                  IF(diff .LT. eps) count=count+1
               END DO
               IF(count .EQ. nc) THEN
                  DO k=1,na
                     l1(k)=beta(bond(k,ia))
                     l2(k)=beta(bond(k,ib))
                  END DO 
                  
*---- Check if label sequence is the same

                  IF(l1(2) .EQ. l2(2) .AND. l1(3) .EQ. l2(3) .OR.
     &                 l1(3) .EQ. l2(2) .AND. l1(2) .EQ. l2(3)) THEN
                     ok=.FALSE.
                     xx(j)=.FALSE.
                     GOTO 200
                  END IF
               END IF
            END DO
200         CONTINUE

*---- Accumulate potentials in work1

            IF(ok) THEN
               nrrbond=nrrbond+1
               work1(nrrbond)=ia
            END IF
         END DO
      END IF      

*---- Write potentials
      
205   DO i=1,nrbond
         ia=work(i)
         IF(nc .EQ. 2 .AND. na .EQ. 2) THEN
            la=bond(1,ia)
            lb=bond(2,ia)
            a=pbond(ia,1)*aux
            b=pbond(ia,2)
            WRITE(kprint,'(2(a7,1x),2f10.3)') beta(la),beta(lb),a,b
         ELSE IF(nc .EQ. 4 .AND. na .EQ. 3) THEN
            la=bond(1,ia)
            lb=bond(2,ia)
            lc=bond(3,ia)
            a=pbond(ia,1)*aux
            b=pbond(ia,2)*aux2
            c=pbond(ia,3)*aux
            d=pbond(ia,4)
            IF(DABS(c) .GT. 1.0D-7 .AND. DABS(d) .GT. 1.0D-7) THEN
               WRITE(kprint,'(3(a7,1x),4f10.3)') beta(la),beta(lb)
     &              ,beta(lc),a,b,c,d
            ELSE
               WRITE(kprint,'(3(a7,1x),2f10.3)') beta(la),beta(lb)
     &              ,beta(lc),a,b
            END IF
         ELSE IF(nc .EQ. 3 .AND. na .EQ. 4) THEN
            la=bond(1,ia)
            lb=bond(2,ia)
            lc=bond(3,ia)
            ld=bond(4,ia)
            IF(pbond(ia,3) .GT. 0) THEN
               a=pbond(ia,1)*aux
               b=pbond(ia,2)*aux2
               WRITE(kprint,'(4(a7,1x),2f10.3,a9)') beta(la)
     &              ,beta(lb),beta(lc),beta(ld),a,b,' HARMONIC'
            ELSE
               a=pbond(ia,1)*aux
               b=pbond(ia,2)
               tmp=IDNINT(pbond(ia,2))
               ntph=IABS(tmp)
               sig=DFLOAT(tmp/ntph)
               IF(a .GT. 0.0D0) c=180.0D0
               IF(a .LE. 0.0D0) c=0.0D0
               a=DABS(a)*sig
               WRITE(kprint,'(4(a7,1x),f10.4,i5,f10.2,a7)')
     &              beta(la),beta(lb),beta(lc),beta(ld),a,ntph,c,
     &              ' COSINE'
            END IF
         END IF
      END DO

      IF(nc .EQ. 2 .AND. na .EQ. 4) THEN
         DO i=1,nrbond
            ia=work(i)
            la=bond(1,ia)
            lb=bond(2,ia)
            lc=bond(3,ia)
            ld=bond(4,ia)
            a=pbond(ia,1)*aux
            b=pbond(ia,2)
            tmp=IDNINT(b)
            ntph=IABS(tmp)
            sig=DFLOAT(tmp/ntph)
            IF(a .GT. 0.0D0) c=180.0D0
            IF(a .LE. 0.0D0) c=0.0D0
            a=DABS(a)*sig
            WRITE(kprint,'(4(a7,1x),f10.4,i5,f10.2)') beta(la)
     &              ,beta(lb),beta(lc),beta(ld),a,ntph,c
         END DO
      END IF


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END


