      SUBROUTINE write_xrms(kprint,n,label,a)

************************************************************************
*   Time-stamp: <97/07/15 18:48:17 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Jun 30 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,kprint
      REAL*8  a(*)
      CHARACTER*2 label


*------------------------- SCRATCH COMMON -----------------------------*

      INTEGER index(m1)
      COMMON /rag1/ index

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,na,nb,map

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      map=0
      DO i=1,n
         IF(DABS(a(i)) .GT. 1.0D-4) THEN
            map=map+1
            index(map)=i
         END IF
      END DO

      IF(MOD(map,6) .EQ. 0) THEN
         DO i=1,map,6
            WRITE(kprint,60000) label,(index(j),a(index(j)),j=i,i+5)
         END DO 
      ELSE
         na=map-MOD(map,6)
         DO i=1,na,6
            WRITE(kprint,60000) label,(index(j),a(index(j)),j=i,i+5)
         END DO
         i=na+1
         nb=MOD(map,6)-1
         WRITE(kprint,60000) label,(index(j),a(index(j)),j=i,i+nb)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

60000 FORMAT(a2,8(1x,i4,1x,f8.3))

      RETURN
      END
