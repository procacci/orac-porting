      SUBROUTINE appgrp(ngrp,grppt,n1,jngrp,jgrppt,shifta,iret,errmsg)

************************************************************************
*                                                                      *
*     This subroutine appends a list of charge group.                  *
*                                                                      *
*                                                                      *
*                                                                      *
*---- Last update 01/10/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley CA 1992                    *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ngrp,jngrp,shifta,n1,iret
      INTEGER grppt(2,*),jgrppt(2,*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m
      LOGICAL err,ok

*==================== EXECUTABLE STATEMENTS ============================

      DO n=1,jngrp
          grppt(1,n+ngrp)=jgrppt(1,n)+shifta
          grppt(2,n+ngrp)=jgrppt(2,n)+shifta
          IF(n1 .LT. ngrp+n) THEN
              errmsg=' Error in APPGRP: Physical dimensions of GRPPT'//
     x               ' exceeded. Abort.'
              iret=1
              RETURN
          END IF
      END DO
      ngrp=ngrp+jngrp

      RETURN
      END
