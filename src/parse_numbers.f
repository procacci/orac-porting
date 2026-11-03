      SUBROUTINE parse_numbers(err_unr,string,nword,list,n_res,nsevere)

************************************************************************
*   Time-stamp: <97/07/08 22:04:39 marchi>                             *
*                                                                      *
*   Read a list of numbers and recognizes token '-': n1 - n2           *
*   meaning that all numbers from n1 to n2 are included in the         *
*   list                                                               *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nword,n_res,nsevere
      CHARACTER*80 string(*)
      INTEGER list(*)
      CHARACTER*80 err_unr(3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER na
      PARAMETER (na=60)
      INTEGER i,naux,naux1,naux2,j
      INTEGER mask(na)
      CHARACTER*80 errmsg
      CHARACTER*8 fmt

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*--- Set mask to 2 for dashes and to 1 otherwise -----------------------
*=======================================================================

      DO i=2,nword
         mask(i)=1
         IF(string(i) .EQ. '-') mask(i)=2
      END DO

*=======================================================================
*---- Try to see if the dash is followed or preceded by another dash ---
*---- In that case abort -----------------------------------------------
*=======================================================================

      DO i=2,nword
         IF(mask(i) .EQ. 2) THEN
            IF(i .EQ. 2 .OR. i .EQ. nword) THEN
               errmsg=err_unr(3)//string(i)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               RETURN
            ELSE IF(mask(i-1) .EQ. 2 .OR. mask(i+1) .EQ. 2) THEN
               errmsg=err_unr(3)//string(i)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
               RETURN
            END IF
         END IF
      END DO

*=======================================================================
*--- mask set 3 for the adjacent strings -------------------------------
*=======================================================================
      
      DO i=2,nword
         IF(mask(i) .EQ. 2) THEN
            mask(i-1)=3
            mask(i+1)=3
         END IF
      END DO

*=======================================================================
*--- Use mask to get the list: different action if mask is 1 or 3 ------
*=======================================================================
      
      naux=0
      DO i=2,nword
         IF(mask(i) .EQ. 1) THEN
            n_res=n_res+1
            CALL fndfmt(1,string(i),fmt)
            READ(string(i),fmt,err=20) list(1+n_res)
         ELSE IF(mask(i) .EQ. 3) THEN
            naux=naux+1
            IF(naux .EQ. 1) THEN
               CALL fndfmt(1,string(i),fmt)
               READ(string(i),fmt,err=20) naux1
            ELSE IF(naux .EQ. 2) THEN
               CALL fndfmt(1,string(i),fmt)
               READ(string(i),fmt,err=20) naux2
               DO j=naux1,naux2
                  n_res=n_res+1
                  list(1+n_res)=j
               END DO
               naux=0
            END IF
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN

 20   CONTINUE
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      END
