      subroutine jartemp(t,temp,tiniz,tfina,velt,rtime,timetiniz
     &     ,timetfina,ucns_h,ucos_h,ucns_l,ucos_l,ucns_m ,ucos_m,urcs_h
     &     ,urcs_l,urcs_m,uslvbon,uslvben,uslvtor ,uslvitor
     &     ,conf_bnd_slv_n1,coul_bnd_slv,coul_bnd_slv_n1 ,self_slv,efact
     &     ,urcp_h,urcp_l,urcp_m,ucop_h,ucop_l ,ucop_m,ucnp_h ,ucnp_l
     &     ,ucnp_m,ubond,ubend,uptors,uitors ,conf_bnd_slt_n1
     &     ,coul_bnd_slt,coul_bnd_slt_n1,self_slt ,ucek,pucek,ucepr,uceh
     &     ,cpress,urcsp_h,urcsp_l,urcsp_m ,ucosp_h,ucosp_l,ucosp_m
     &     ,ucnsp_h,ucnsp_l,ucnsp_m,eer_m,eer_l,eer_h,volume,pext,gascon
     &     ,ninn0,nplottempstn0,tn0,kplot_steer_temp,nsteer_temp)

************************************************************************
*                                                                      *
*     JARTEMP  irreversibly changes the external bath temperature and  *
*     compute the adimensional work in the noneq transformation        *
*     according to the Generalized Crooks theorem by Chelly et al.     *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, Florence Univ 2008                    *
*                                                                      *
*______________________________________________________________________*



*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 t,temp,tiniz,tfina,velt,rtime,timetiniz,timetfina,ucns_h
     &     ,ucos_h,ucns_l,ucos_l,ucns_m ,ucos_m,urcs_h ,urcs_l,urcs_m
     &     ,uslvbon,uslvben,uslvtor ,uslvitor ,conf_bnd_slv_n1
     &     ,coul_bnd_slv,coul_bnd_slv_n1 ,self_slv,efact,urcp_h,urcp_l
     &     ,urcp_m,ucop_h,ucop_l ,ucop_m,ucnp_h ,ucnp_l,ucnp_m,ubond
     &     ,ubend,uptors,uitors ,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt ,ucek,pucek,ucepr,uceh,urcsp_h
     &     ,urcsp_l,urcsp_m,ucosp_h,ucosp_l,ucosp_m,ucnsp_h,ucnsp_l
     &     ,ucnsp_m ,eer_m,eer_l,eer_h,volume,pext,gascon,tn0

      INTEGER ninn0,nplottempstn0,kplot_steer_temp,nsteer_temp

      LOGICAL cpress
      
*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 ektot,hextended,worktemp,ktbeta,ustot,uptot,upstot,dt
      SAVE   hextended,worktemp

*==================== EXECUTABLE STATEMENTS ============================

c      write(6,*) t,temp,volume,tn0,nplottempstn0,kplot_steer_temp 
c      write(6,*) timetiniz,timetfina,rtime
      
      dt=abs(rtime-timetfina)
      if(rtime.ge.timetiniz.and.rtime.lt.timetfina) THEN  
        nsteer_temp=nsteer_temp+1
        if(nsteer_temp.eq.1) THEN 
          write(kplot_steer_temp,1000) 
1000      FORMAT('# Adimensional work using the Generalized '
     &         'Crooks theorem for non equilibrium thermal changes'//,
     &         '#',7x,'time/fs',4x,'Text/K',4x,'Tint/K',9x,
     &         'Omega',3x,'Hext/KJmol-1')
        end if
        t = tiniz + velt*(rtime-timetiniz)
        ustot=(ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m+urcs_h
     &       +urcs_l+urcs_m+uslvbon+uslvben+uslvtor+uslvitor
     &       +conf_bnd_slv_n1+coul_bnd_slv+coul_bnd_slv_n1
     &       +self_slv)*efact/1000.0D0
        uptot=(urcp_h+urcp_l+urcp_m+ucop_h+ucop_l+ucop_m+ucnp_h
     &       +ucnp_l+ucnp_m+ubond+ubend+uptors+uitors
     &     +conf_bnd_slt_n1+coul_bnd_slt+coul_bnd_slt_n1
     &       +self_slt)*efact/1000.d0
        
        ektot = (ucek+pucek)/1000.d0  
      
        IF(cpress) ektot=ektot+ucepr/1000.0D0
        ektot=ektot+uceh/1000.0D0
        
        upstot= 0.d0
        upstot = (urcsp_h+urcsp_l+urcsp_m+ucosp_h+ucosp_l+ucosp_m
     &       +ucnsp_h+ucnsp_l+ucnsp_m)*efact/1000.d0
        
c--           According to the GCT The extended system energy DOES NOT 
c--           include the hpot POTENTIAL bath energy   
c       write(6,*) ektot,upstot,ustot,uptot
        Hextended= ektot+ustot+uptot+upstot+(eer_m+eer_l+eer_h)
     &       *efact/1000.d0
        if(cpress) THEN 
          Hextended=Hextended+(pext*volume*efact)/1000.0D0
        end if
        
        Hextended=ustot+uptot+upstot+ektot
c        write(6,*) "gascon",gascon
        ktbeta=1/(t*gascon/1000.d0) 
c--           No dimession for worktemp

        worktemp=worktemp-ktbeta*velt*Hextended*tn0/t 
c        write(6,*) "gascon",gascon,nplottempstn0
        if(mod(ninn0,nplottempstn0).EQ.0.or.dt.LE.tn0) THEN 
          write(kplot_steer_temp,1015) rtime,t,temp,worktemp,Hextended 
1015      FORMAT(f15.1,2f10.2,2f15.4)
        end if 
      end if
      return
      end
