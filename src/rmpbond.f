      SUBROUTINE rmpbond(stretch,stretch_heavy,lbnd,lbond,lstrtch
     &     ,lstretch,lcnstr,lconstr,n1,type,pbondx,pbondy)

************************************************************************
*   Time-stamp: <97/12/10 15:28:18 marchi>                             *
*                                                                      *
*   Remove bond stretching and constraint when necessary               *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Jun 23 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER      lbond,lstretch,lconstr,n1,lbnd(2,*),lstrtch(2,*)
     &     ,lcnstr(2,*)
      REAL*8       pbondx(n1,*),pbondy(n1,*)
      CHARACTER*7  type(*)
      LOGICAL stretch_heavy,stretch
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER      i,ncount1,ncount2,i1,i2

*==== EXECUTABLE STATEMENTS: ==========================================*

      IF(stretch) THEN
         lstretch=lbond
         lconstr=0
         DO i=1,lbond
            lstrtch(1,i)=lbnd(1,i)
            lstrtch(2,i)=lbnd(2,i)
         END DO
         IF(stretch_heavy) THEN
            ncount2=0
            DO i=1,lstretch
              i1=1
              i2=1
              if(type(lstrtch(1,i))(1:1).eq."*") i1=2
              if(type(lstrtch(2,i))(1:1).eq."*") i2=2
               IF(  type(lstrtch(1,i))(i1:i1) .EQ. 'h' .OR.
     &              type(lstrtch(2,i))(i2:i2) .EQ. 'h' ) THEN
                  ncount2=ncount2+1
                  lcnstr(1,ncount2)=lstrtch(1,i)
                  lcnstr(2,ncount2)=lstrtch(2,i)
                  IF(DABS(pbondx(i,1)) .LT. 1.0D-7) THEN
                     pbondy(ncount2,1)=500.0D0
                  ELSE 
                     pbondy(ncount2,1)=pbondx(i,1)
                  END IF
                  pbondy(ncount2,2)=pbondx(i,2)
               END IF
            END DO
            lconstr=ncount2

            ncount1=0
            DO i=1,lstretch
              i1=1
              i2=1
              if(type(lstrtch(1,i))(1:1).eq."*") i1=2
              if(type(lstrtch(2,i))(1:1).eq."*") i2=2
               IF(  type(lstrtch(1,i))(i1:i1) .NE. 'h' .AND.
     &              type(lstrtch(2,i))(i2:i2) .NE. 'h' ) THEN
                  ncount1=ncount1+1
                  lstrtch(1,ncount1)=lstrtch(1,i)
                  lstrtch(2,ncount1)=lstrtch(2,i)
                  pbondx(ncount1,1)=pbondx(i,1)
                  pbondx(ncount1,2)=pbondx(i,2)
               END IF
            END DO
            lstretch=ncount1
         END IF
      ELSE
         lstretch=0
         lconstr=lbond
         DO i=1,lbond
            lcnstr(1,i)=lbnd(1,i)
            lcnstr(2,i)=lbnd(2,i)
            pbondy(i,1)=1.0D0
            pbondy(i,2)=pbondx(i,2)
         END DO
      END IF
         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
