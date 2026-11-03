      FUNCTION f1dim_total(mapnl,mapdn,nmapdn
     &     ,nnlpp0,nnlpp,fudgec,x,xp0,yp0,zp0,fpx,fpy,fpz)

************************************************************************
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      use parst 
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER mapnl(*),mapdn(2,*),nmapdn(*)
      INTEGER nnlpp0(*),nnlpp(*)
      REAL*8  fudgec,f1dim_total
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),x

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i
      REAL*8  utotal

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------
  

*==================== EXECUTABLE STATEMENTS ============================


      DO i=1,ntap
        if(mass(i).lt.massinfty) THEN 
          xp0(i)=xp0(i)+x*fpx(i)
          yp0(i)=yp0(i)+x*fpy(i)
          zp0(i)=zp0(i)+x*fpz(i)
        endif
      END DO
      CALL comp_total_energy(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,x,xp0,yp0,zp0,utotal)
#ifndef _PGI_
      if(isnan(utotal)) THEN
        write(6,*) "NaN utotal in f1dim_total; progran stops"
        STOP
      ENDIF
#endif
      DO i=1,ntap
        if(mass(i).lt.massinfty) THEN 
         xp0(i)=xp0(i)-x*fpx(i)
         yp0(i)=yp0(i)-x*fpy(i)
         zp0(i)=zp0(i)-x*fpz(i)
       endif
       END DO
      f1dim_total=utotal
      RETURN
      END

