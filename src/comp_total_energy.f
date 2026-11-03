      SUBROUTINE comp_total_energy(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,x,xp0,yp0,zp0,utotal)

************************************************************************
*   Time-stamp: <98/03/06 17:53:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 14 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER mapnl(*),mapdn(2,*),nmapdn(*)
      INTEGER nnlpp0(*),nnlpp(*)
      REAL*8  fudgec

      REAL*8  xp0(*),yp0(*),zp0(*),x,prt(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors,uconf,ucoul
     &     ,ureal,urecp,urcs,urcp,urcsp,urcsp_m,eer,ucns,ucos,ucnsp
     &     ,ucosp,ucnp,ucop,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,virs,virsp,virp,gsin14,fsbend,fsbond,gsbend
     &     ,gsbond,unb14,cnb14,fsin14,ungrp,cngrp,uumb,gr,ubond,uslvbon
     &     ,uslvben,unbond,cnbond,purecp,upconf,upcoul,enisolv,press
     &     ,stressd(3,3),stressr(3,3),st(3,3),pv,ustot,uptot,upstot
     &     ,utotal,fpx,fpy,fpz
      INTEGER iret,abmd_dir,nstep
      CHARACTER*80 errmsg

      CHARACTER*1 rshell,rshk
      DATA   rshell/'h'/rshk/'h'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      nstep=0
      CALL get_total_energy(.FALSE.,mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,abmd_dir,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,prt,utotal
     &     ,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,fsin14,unb14,cnb14,uslvbon,uslvben,uslvtor,uslvitor
     &     ,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer,enisolv,nstep)

        
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
