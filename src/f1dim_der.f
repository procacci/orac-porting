      FUNCTION f1dim_der(iatom,icoord,mapnl,mapdn
     &     ,nmapdn,nnlpp0,nnlpp,fudgec,h,xp0,yp0,zp0,fpx,fpy
     &     ,fpz)

************************************************************************
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use parst
      use cpropar
      use unit, only:efact

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iatom,icoord,mapnl(*),mapdn(2,*),nmapdn(*)
      INTEGER nnlpp0(*),nnlpp(*)
      REAL*8  fudgec,f1dim_der
      REAL*8  h,xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)


*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,nstep
      REAL*8  prt(3,3),utotal,ucns,ucos,urcs,coul_bnd_slv
     &     ,conf_bnd_slv_n1,coul_bnd_slv_n1,self_slv,uslvbon,uslvben
     &     ,uslvtor,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop
     &     ,urcp,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt
     &     ,ucnsp,ucosp,urcsp,eer,fsin14,unb14,cnb14,enisolv

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------
  

*==================== EXECUTABLE STATEMENTS ============================

      nstep=0
      IF(icoord .EQ. 1) THEN
         xp0(iatom)=xp0(iatom)+h
      ELSE IF(icoord .EQ. 2) THEN
         yp0(iatom)=yp0(iatom)+h
      ELSE IF(icoord .EQ. 3) THEN
         zp0(iatom)=zp0(iatom)+h
      END IF

      CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,1,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,prt,utotal
     &     ,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,fsin14,unb14,cnb14,uslvbon,uslvben,uslvtor,uslvitor
     &     ,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer,enisolv,nstep)

      IF(icoord .EQ. 1) THEN
         xp0(iatom)=xp0(iatom)-h
      ELSE IF(icoord .EQ. 2) THEN
         yp0(iatom)=yp0(iatom)-h
      ELSE IF(icoord .EQ. 3) THEN
         zp0(iatom)=zp0(iatom)-h
      END IF
      if(agbnp.OR.lsasa) THEN 
        utotal=utotal+enisolv*efact/1000.d0
      END IF
      f1dim_der=utotal
      
      RETURN
      END

