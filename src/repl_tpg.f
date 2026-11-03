      SUBROUTINE repl_tpg(mapnl,gnmol,iret,errmsg)

************************************************************************
*   Time-stamp: <98/03/19 17:04:59 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Mar 14 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER gnmol,iret,mapnl(*)
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER idummy,i

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL dupl_int2(.TRUE.,ntap,lacc,llacc,2,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for lacc. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,ldon,lldon,2,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for ldon. Abort'
         RETURN
      END IF

      CALL dupl_int1(.FALSE.,idummy,nbtype,ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for nbtype. Abort'
         RETURN
      END IF

      CALL dupl_int1(.TRUE.,nbun,nres(1,1),ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for nres. Abort'
         RETURN
      END IF

      CALL dupl_int1(.FALSE.,idummy,nres(1,2),ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for nres. Abort'
         RETURN
      END IF

      CALL dupl_int1(.TRUE.,ntap,mback,ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for mback. Abort'
         RETURN
      END IF

      CALL dupl_int1(.FALSE.,idummy,mend,nbun,nores,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for mend. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,grppt,ngrp,2,m11,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for grppt. Abort'
         RETURN
      END IF
      CALL dupl_int2(.TRUE.,ntap,lbnd,lbond,2,m9,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for lbnd. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,lbndg,lbend,3,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for lbndg. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,ltor,ltors,4,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for ltor. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,litr,litor,4,m4,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for litr. Abort'
         RETURN
      END IF
 
      CALL dupl_int2(.TRUE.,ntap,int14,int14p,2,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for int14. Abort'
         RETURN
      END IF

      CALL dupl_int1(.FALSE.,ntap,type14,int14p,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for type14. Abort'
         RETURN
      END IF

      CALL dupl_int2(.TRUE.,ntap,int13,int13p,2,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for int13. Abort'
         RETURN
      END IF

      CALL dupl_concta(ntap,concta,ntap,m10,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for concta. Abort'
         RETURN
      END IF

      CALL dupl_mapnl(ntap,mapnl,ntap,m8,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for mapnl. Abort'
         RETURN
      END IF

      CALL dupl_real1(chrge,ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for chrge. Abort'
         RETURN
      END IF

      CALL dupl_real1(mass,ntap,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for mass. Abort'
         RETURN
      END IF

****  Need to do arrays type potbe(i,2) e concta(i,j) type of arrays and 
****  nrigg, prsymb, mend and character arrays

      CALL dupl_reali2(potbo,lbond,2,m9,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for potbo. Abort'
         RETURN
      END IF

      CALL dupl_reali2(potbe,lbend,4,m2,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for potbe. Abort'
         RETURN
      END IF

      IF(amphi) THEN
         CALL dupl_reali2(ptorj,ltors,4,m2,gnmol,iret)
         IF(iret .EQ. 1) THEN
            errmsg='In REPL_TPG: Dimension exceeded for ptorj. Abort'
            RETURN
         END IF
      ELSE
         CALL dupl_reali2(potto,ltors,2,m2,gnmol,iret)
         IF(iret .EQ. 1) THEN
            errmsg='In REPL_TPG: Dimension exceeded for ptorj. Abort'
            RETURN
         END IF
      END IF
      
      CALL dupl_reali2(potit,litor,3,m4,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for ptorj. Abort'
         RETURN
      END IF

      IF(nrigg .NE. 0) THEN
         iret=1
         errmsg='In REPL_TPG: Rigid atomic groups are incompatible '
     &//' with symmetry. Abort'
         RETURN
      END IF
      
      CALL dupl_char1(beta,ntap,7,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for beta. Abort'
         RETURN
      END IF

      CALL dupl_char1(betb,ntap,7,m1,gnmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In REPL_TPG: Dimension exceeded for betb. Abort'
         RETURN
      END IF


      nbun=nbun*gnmol
      llacc=llacc*gnmol
      lldon=lldon*gnmol
      nbone=nbone*gnmol
      ntap=ntap*gnmol
      lbond=lbond*gnmol
      lbend=lbend*gnmol
      ltors=ltors*gnmol
      litor=litor*gnmol
      int14p=int14p*gnmol
      int13p=int13p*gnmol
      ngrp=ngrp*gnmol
      wmtp=wmtp*DFLOAT(gnmol)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
