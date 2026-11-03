      SUBROUTINE dumptp(kdump,mapnl,nmapnl)

************************************************************************
*                                                                      *
*                                                                      *
*     DUMPRS will dump all data necessary to produce a restart         *
*     to the data file.                                                *
*                                                                      *
*                                                                      *
*----- Last update 06/5/90 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS none.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER kdump,nmapnl,mapnl(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j

*==================== EXECUTABLE STATEMENTS ============================

      REWIND kdump

      WRITE(kdump,*) llacc,lldon,nbone,ntap,lbond,lbend
     &     ,ltors,litor,lpnbd,int14p,int13p,nspw,nunit,lphyd,ngrp,nrigg
     &     ,nrigat,nmapnl,nmol,nato_slv,adihed,slv_exist,slt_exist
      WRITE(kdump,*) ss_point(1,1),ss_point(1,2)
      WRITE(kdump,*) (ss_point(i+1,1),i=1,ss_point(1,1)),(ss_point(i+1
     &     ,2),i=1,ss_point(1,2))
      WRITE(kdump,*) (ss_index(i),i=1,ntap)
      WRITE(kdump,*) (lacc(1,i),lacc(2,i),i=1,llacc),
     x             (ldon(1,i),ldon(2,i),i=1,lldon)
      WRITE(kdump,*) (nbtype(i),i=1,ntap)
      WRITE(kdump,*) (nres(i,1),i=1,ntap)
      WRITE(kdump,*) (nres(i,2),i=1,ntap) 
      WRITE(kdump,*) (mback(i),i=1,ntap) 
      WRITE(kdump,*) (grppt(1,i),i=1,ngrp),(grppt(2,i),i=1,ngrp)
      WRITE(kdump,*) (lbnd(1,i),i=1,lbond),(lbnd(2,i),i=1,lbond)
      WRITE(kdump,*)     (concta(i,1),(concta(i,j),j=2,concta(i,1)+1),i
     &     =1,ntap)
      WRITE(kdump,*) ((lbndg(j,i),i=1,lbend),j=1,3)
      WRITE(kdump,*) ((ltor(j,i),i=1,ltors),j=1,4)
      WRITE(kdump,*) ((litr(j,i),i=1,litor),j=1,4)
      WRITE(kdump,*) ((int14(j,i),j=1,2),i=1,int14p),
     x             (type14(i),i=1,int14p),
     x             ((int13(j,i),j=1,2),i=1,int13p),
     x             ((nhtype(i,j),j=1,lpnbd),i=1,lpnbd)
      WRITE(kdump,*) wmtp,(chrge(i),i=1,ntap),(mass(i),i=1,ntap)
      WRITE(kdump,*) (mapnl(i),i=1,nmapnl)
      IF(amphi) THEN
          WRITE(kdump,*) nbun,iz
          WRITE(kdump,*) ((potbo(i,j),i=1,lbond),j=1,2),
     x                 ((potbe(i,j),i=1,lbend),j=1,7),
     x                 ((ptorj(i,j),i=1,ltors),j=1,4),
     x                 ((potit(i,j),i=1,litor),j=1,3),
     x                 (pnbd1(i),i=1,lpnbd),(pnbd2(i),i=1,lpnbd),
     x                 (pnbd3(i),i=1,lpnbd),(pnbd4(i),i=1,lpnbd),
     &                 (sjorg(i),i=1,lpnbd),
     x                 (ejorg(i),i=1,lpnbd),(phyd1(i),i=1,lphyd),
     x                 (phyd2(i),i=1,lphyd)
          IF(iz .EQ. 1) THEN
             WRITE(kdump,*) 
     &                 (ecc6(i),i=1,lpnbd*(lpnbd+1)/2),
     &                 (ecc12(i),i=1,lpnbd*(lpnbd+1)/2)
          END IF
      ELSE
          WRITE(kdump,*) nbun,iz
          WRITE(kdump,*) ((potbo(i,j),i=1,lbond),j=1,2),
     x                 ((potbe(i,j),i=1,lbend),j=1,7),
     x                 ((potto(i,j),i=1,ltors),j=1,2),
     x                 ((potit(i,j),i=1,litor),j=1,3),
     x                 (pnbd1(i),i=1,lpnbd),(pnbd2(i),i=1,lpnbd),
     x                 (pnbd3(i),i=1,lpnbd),(pnbd4(i),i=1,lpnbd),
     &                 (phyd1(i),i=1,lphyd),
     x                 (phyd2(i),i=1,lphyd)
          IF(iz .EQ. 1) THEN
             WRITE(kdump,*) 
     &                 (ecc6(i),i=1,lpnbd*(lpnbd+1)/2),
     &                 (ecc12(i),i=1,lpnbd*(lpnbd+1)/2)
          END IF
      END IF
      WRITE(kdump,*) (nasitp(i),i=1,nrigg), (nbsitp(i),i=1,nrigg)
      WRITE(kdump,*) ((asitp(j,i),j=1,nasitp(i)),i=1,nrigg)
      WRITE(kdump,*) ((bsitp(j,i),j=1,nbsitp(i)),i=1,nrigg)
      WRITE(kdump,'(12A7)') (beta(i),i=1,ntap),(betb(i),i=1,ntap),
     x             (alnbd(i),i=1,lpnbd)
      WRITE(kdump,'(12A8)') (prsymb(i),i=1,nunit)
      WRITE(kdump,*) (mend(i),i=1,nbun)
      CLOSE(kdump)

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
