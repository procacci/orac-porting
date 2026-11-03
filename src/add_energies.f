      SUBROUTINE add_energies(pme,pressure,slv_exist,slt_exist,ucns,ucos
     &     ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1,self_slv
     &     ,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors,uitors,ubond
     &     ,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp,eer,uconf,ucoul
     &     ,ureal,urecp,pubnd,purecp,puconf,pucoul,upconf,upcoul,stressd
     &     ,stressr,stress_conf,stress_tot,co,oc,volume,unitp,press_conf
     &     )

************************************************************************
*   Time-stamp: <95/01/07 00:47:53 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 22 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL pme,slv_exist,slt_exist,pressure
      REAL*8  ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,self_slv,uslvbon,uslvben,uslvtor,uslvitor
     &     ,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer,uconf,ucoul,ureal,urecp,pubnd,purecp,puconf
     &     ,pucoul,upconf,upcoul,stressd,stressr,stress_conf,stress_tot
     &     ,co,oc,volume,unitp,press_conf

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  unbond,cnbond

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(slv_exist) THEN
         uconf=ucns+ucos+urcs+coul_bnd_slv+conf_bnd_slv_n1
     &        +coul_bnd_slv_n1+self_slv + uumb
         ucoul=ucos+urcs+coul_bnd_slv+coul_bnd_slv_n1+self_slv
         ureal=ucos+coul_bnd_slv_n1
         urecp=urcs+self_slv+coul_bnd_slv
      END IF
      
*---        Energy terms for the protein         -----------------------
      
      IF(slt_exist) THEN 
         pubnd=uptors+uitors+ubond+ubend
         unbond=ucnp
         cnbond=ucop
         purecp=urcp
         puconf= unbond + conf_bnd_slt_n1 + uumb
         pucoul= cnbond + purecp + coul_bnd_slt + coul_bnd_slt_n1
     &        +self_slt
      END IF   
      
*---        Mixed terms                          -----------------------
      
      IF(slv_exist .AND. slt_exist) THEN 
         upconf= ucnsp + ucosp + urcsp
         upcoul= ucosp + urcsp
      END IF
      
      IF(pressure) THEN
         CALL comp_stress_conf(stressd,stressr,stress_conf,oc
     &        ,volume,unitp,press_conf)
         CALL dcopy(9,stress_conf,1,stress_tot,1)
      END IF
      
*----       If pme=T no solute-solvent separation ----------------------
      
      IF(pme) THEN 
         IF(slv_exist.AND.(.NOT.slt_exist)) THEN
            uconf=uconf + eer
            ucoul=ucoul + eer
            urecp= eer 
         END IF
         IF(slt_exist.AND.(.NOT.slv_exist)) THEN
            purecp=eer 
            pucoul=pucoul+eer 
         END IF
         IF(slt_exist.AND.slv_exist) THEN 
            upconf=upconf + eer 
            upcoul=upcoul + eer 
         END IF
      END IF
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      
      RETURN
      END
