      SUBROUTINE add_solvent_tpg(mapnl,mapnl_slv,iret,errmsg)

************************************************************************
*   Time-stamp: <98/03/17 13:40:38 marchi>                             *
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
      use unit, only:kprint
      
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,mapnl(*),mapnl_slv(*),ndimr,j,nsevere
      CHARACTER*132 errmsg,string

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER idummy,i
      REAL*8  sum

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      nsevere=0
      CALL add_int2(.TRUE.,ntap,nato_slv,lacc,lacc_slv,llacc,llacc_slv,2
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=llacc+llacc_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)'Dimension exceeded for lacc.' //
     &        ' increase _SIT_SOLU_ in config.H to ' // string(1:j)
         nsevere=nsevere+1
101      FORMAT(A150)
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,ldon,ldon_slv,lldon,lldon_slv,2
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=lldon+lldon_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for ldon:'//
     &        ' increase _SIT_SOLU_ in config.H to ' //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int1(.FALSE.,ntap,nato_slv,ntap,nato_slv,nbtype
     &     ,nbtype_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
        errmsg
     &     =' Dimension exceeded for nbtype:'//
     &        ' increase _SIT_SOLU_ in config.H to ' //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int1(.TRUE.,nbun,nbun_slv,ntap,nato_slv,nres(1,1)
     &     ,nres_slv(1,1),m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for nres:'//
     &        ' increase _SIT_SOLU_ in config.H to ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_int1(.FALSE.,nbun,nbun_slv,ntap,nato_slv,nres(1,2)
     &     ,nres_slv(1,2),m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for nres:'//
     &        ' increase _SIT_SOLU_ in config.H to '  //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int1(.TRUE.,ntap,nato_slv,ntap,nato_slv,mback,mback_slv
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap*nato_slv
        call int_str(ndimr,string,j)
         errmsg
     &        =' Dimension exceeded for mback:'//
     &        ' increase _SIT_SOLU_ in config.H to '  //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int1(.FALSE.,idummy,idummy,nbun,nbun_slv,mend,mend_slv
     &     ,nores,nmol,iret)
      IF(iret .EQ. 1) THEN
         write(kprint,101)' Dimension exceeded for NRES. Abort'
         write(kprint,102) nores,nbun,nbun_slv,nmol, nbun+nbun_slv*nmol
102      FORMAT(" == _NRES_ in config.H is seto to ",i5 / 
     &          " == solute units",I5 /
     &          " == residues defining solvent",I5 /
     &          " == solvent molecules",I5 /
     &          " Action: increase _NRES_ in config.H to at least", I5 )
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,grppt,grppt_slv,ngrp,ngrp_slv
     &     ,2,m11,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ngrp+ngrp_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &        =' Dimension exceeded for grptt:'//
     &        ' increase _TGROUP_ in config.H to ' // string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,lbnd,lbnd_slv,lbond,lbond_slv,2
     &     ,m9,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=lbond+lbond_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for lbnd:'//
     &        ' increase _SIT_SOLU_ in config.H to ' // string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,lbndg,lbndg_slv,lbend,lbend_slv
     &     ,3,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=lbend+lbend_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for lbndg:'//
     &        ' increase _SIT_SOLU_ in config.H to '  // string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,ltor,ltor_slv,ltors,ltors_slv,4
     &     ,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ltors+ltors_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for ltor:'//
     &   ' increase _SIT_SOLU_ in config.H or change m2 in parst.f90'//
     &   ' Torsions are ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,litr,litr_slv,litor,litor_slv,4
     &     ,m4,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=litor+litor_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for litr:'//
     &  ' increase _SIT_SOLU_ in config.H or change m4 in parst.f90'  //
     &  ' I-Torsions are ' //string(1:j) 
         nsevere=nsevere+1
      END IF
 
      CALL add_int2(.TRUE.,ntap,nato_slv,int14,int14_slv,int14p
     &     ,int14p_slv,2,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=int14p+int14p_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for int14:'//
     &   ' increase _SIT_SOLU_ in config.H or change m2 in parst.f90' //
     &   ' int14 are ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_int1(.FALSE.,ntap,nato_slv,int14p,int14p_slv,type14
     &     ,type14_slv,m2,nmol,iret)

      IF(iret .EQ. 1) THEN
        ndimr=int14p+int14p_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for type14:'//
     &   ' increase _SIT_SOLU_ in config.H or change m2 in parst.f90' //
     &   ' int14p are ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,int13,int13_slv,int13p
     &     ,int13p_slv,2,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=int13p+int13p_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for int13:'//
     &   ' increase _SIT_SOLU_ in config.H or change m2 in parst.f90' //
     &   ' int13 are ' //string(1:j) 
         nsevere=nsevere+1
      END IF

*****Stopped here for the moment ******************************

      CALL add_concta(ntap,nato_slv,concta,concta_slv,ntap,nato_slv,m10
     &     ,m1,slvatm,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for concta:' //
     &       ' increase _SIT_SOLU_ in config.H to ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_mapnl(ntap,nato_slv,ntap,nato_slv,mapnl,mapnl_slv,m8,nmol
     &     ,iret)
      IF(iret .EQ. 1) THEN
        call int_str(m8,string,j)
         errmsg
     &    =' Dimension exceeded for mapnl:' //
     &  ' increase _SIT_SOLU_ in config.H or change m8 in parst.f90.' //
     &  ' m8 is ' //string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_real1(ntap,nato_slv,chrge,chrge_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for chrge:' //
     &   ' increase _SIT_SOLU_ in config.H to ' //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_real1(ntap,nato_slv,mass,mass_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for mass:'//
     &   ' increase _SIT_SOLU_ in config.H to ' // string(1:j)
         nsevere=nsevere+1
      END IF

****  Need to do arrays type potbe(i,2) e concta(i,j) type of arrays and 
****  nrigg, prsymb, mend and character arrays

      CALL add_real2(2,lbond,lbond_slv,potbo,potbo_slv,m9,slv2,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
        ndimr=lbond+lbond_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &   =' Dimension exceeded for potbo:'//
     &   ' increase _SIT_SOLU_ in config.H or change m9 in parst.f90' //
     &   ' req dim is '  //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_real2(7,lbend,lbend_slv,potbe,potbe_slv,m2,slv3,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
        ndimr=lbend+lbend_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &   =' Dimension exceeded for potbe:'//
     &   ' increase _SIT_SOLU_ in config.H or change m2 in parst.f90.'//
     &   ' req. dim is '  //string(1:j)
         nsevere=nsevere+1
      END IF

      CALL add_real2(2,ltors,ltors_slv,potto,potto_slv,m3,slv4,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
        ndimr=ltors+ltors_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for ptorj:'//
     &   ' increase _SIT_SOLU_ in config.H or change m3 in parst.f90.'//
     &   ' req. dim is '  //string(1:j)
         nsevere=nsevere+1
      END IF
      
      CALL add_real2(3,litor,litor_slv,potit,potit_slv,m4,slv5,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
        ndimr=litor+litor_slv*nmol
        call int_str(ndimr,string,j)
         errmsg
     &    =' Dimension exceeded for ptorj:'//
     &   ' increase _SIT_SOLU_ in config.H or change m4 in parst.f90.'//
     &   ' m4 is '  //string(1:j) 
      END IF

      IF(nrigg .NE. 0) THEN
         iret=1
         errmsg
     &        =' Rigid atomic groups are '
     &        //'incompatible with symmetry. Abort'
         nsevere=nsevere+1
      END IF
      
      CALL add_char1(7,ntap,nato_slv,beta,beta_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for beta:'//
     &   ' increase _SIT_SOLU_ in config.H to '//string(1:j) 
         nsevere=nsevere+1
      END IF

      CALL add_char1(7,ntap,nato_slv,betb,betb_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
        ndimr=ntap+nato_slv*nmol
        call int_str(ndimr,string,j)
         write(kprint,101)' Dimension exceeded for betb:'//
     &   ' increase _SIT_SOLU_ in config.H to '// string(1:j) 
         nsevere=nsevere+1
      END IF
      if(nsevere.gt.0) return

      nbun=nbun+nbun_slv*nmol
      llacc=llacc+llacc_slv*nmol
      lldon=lldon+lldon_slv*nmol
      nbone=nbone+nbone_slv*nmol
      ntap=ntap+nato_slv*nmol
      lbond=lbond+lbond_slv*nmol
      lbend=lbend+lbend_slv*nmol
      ltors=ltors+ltors_slv*nmol
      litor=litor+litor_slv*nmol
      int14p=int14p+int14p_slv*nmol
      int13p=int13p+int13p_slv*nmol
      ngrp=ngrp+ngrp_slv*nmol
      sum=0.0D0
      DO i=1,nato_slv
         sum=sum+mass_slv(i)
      END DO
      wmtp=wmtp+DFLOAT(nmol)*sum

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
