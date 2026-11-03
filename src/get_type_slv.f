      SUBROUTINE get_type_slv(nato_slt,nato_slv,beta,beta_slv,nbeta_slv
     &     ,ntype_slv,type_slv,itypes,offset,types_gofr,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/02 14:36:12 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Apr 20 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato_slv,nato_slt,ntype_slv,type_slv(*),itypes(*),offset
     &     ,nbeta_slv,types_gofr,iret
      CHARACTER*7 beta(*)
      CHARACTER*1 beta_slv(*)
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,iox,count

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(nato_slv .EQ. 0) THEN
         ntype_slv=0
         offset=0
         RETURN 
      END IF
      
      DO i=1,nbeta_slv
         beta_slv(i)=' '
         itypes(i)=0
      END DO
      count=0
      DO i=nato_slt+1,nato_slt+nato_slv
         iox=0
         IF(count .NE. 0) THEN
            DO j=1,count
               IF(beta(i)(1:1) .EQ. beta_slv(j)) THEN
                  iox=1
                  itypes(j)=itypes(j)+1
                  type_slv(i-nato_slt)=j
               END IF
            END DO
         ENDIF 
         IF(iox .EQ. 0) THEN
            count=count+1
            IF(types_gofr .LT. count) THEN
               iret=1
               errmsg=
     &' Size of the gofrs arrays too small. Increase _TYP_SOLV_GOFR_'
     &//' in config.H.'
               RETURN
            END IF
            beta_slv(count)=beta(i)(1:1)
            itypes(count)=itypes(count)+1
            type_slv(i-nato_slt)=count
         END IF
      END DO
      ntype_slv=count
      offset=ntype_slv*(ntype_slv+1)/2
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
