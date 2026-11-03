      SUBROUTINE redefine_slv(ss_point,m1,ss_index,prsymb,mend,nbun
     &     ,slv_redef_type,mres,grppt,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/15 20:07:10 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Jun 27 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER m1,iret
      INTEGER nbun,ss_point(m1,*),ss_index(*),mend(*),mres(2,*),grppt(2,
     &     *)
      CHARACTER*8 prsymb(*),slv_redef_type
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ia,ib,i1,map1,map2,m
      LOGICAL ok,something_wrong

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      something_wrong=.TRUE.
      DO ia=1,nbun
         IF(prsymb(mend(ia))(1:8) .EQ. slv_redef_type(1:8)) THEN
            ok=.TRUE.
            something_wrong=.FALSE.
         ELSE
            ok=.FALSE.
         END IF
         IF(ok) THEN
            DO ib=mres(1,ia),mres(2,ia)
               DO i1=grppt(1,ib),grppt(2,ib)
                  ss_index(i1)=2
               END DO
            END DO
         END IF
      END DO
      IF(something_wrong) THEN
         iret=2
         errmsg=' Solvent type '//slv_redef_type//
     &        ' not found. Solvent not redefined.'
         RETURN
      END IF
         
      map1=0
      map2=0
      DO ia=1,nbun
         DO ib=mres(1,ia),mres(2,ia)
            DO i1=grppt(1,ib),grppt(2,ib)
               IF(ss_index(i1) .EQ. 1) THEN
                  map1=map1+1
                  ss_point(map1+1,1)=i1
               ELSE IF(ss_index(i1) .EQ. 2) THEN
                  map2=map2+1
                  ss_point(map2+1,2)=i1
               END IF
            END DO
         END DO
      END DO
      ss_point(1,1)=map1
      ss_point(1,2)=map2

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
