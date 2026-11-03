      SUBROUTINE bldbox(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,iret
     &     ,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),gh(*)
      CHARACTER*80 errmsg
      INTEGER mapnl(*),mapnl_slv(*)


*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nb,nc,i
      PARAMETER (nb=m1*5,nc=192)
      INTEGER list(3,nb),hlist(3,nb),rlist(3,nb),index(nb)
      REAL*8  xpa(nb),ypa(nb),zpa(nb),xpb(nb),ypb(nb),zpb(nb),cg(3,3,nc)
     &     ,tg(3,nc),work(3*m1)
      COMMON / rag1 / xpa,ypa,zpa,xpb,ypb,zpb,cg,tg,work,list,hlist
     &     ,rlist,index
      INTEGER m,ib,i1
      
      IF(nb .LT. (nmol+1)*nato_slv) THEN
        write(6,'(A3,I6,A6,I6,A9,I6)') "NB=", nb, " NMOL=",nmol
     &       ," NATO_SLV=", nato_slv
        errmsg=' IN BLDBOX : Insufficient dimensions do not allow to'
     &//' generate solvent.'
          iret=1
          RETURN
      END IF

      CALL setup_box(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,xpa,ypa
     &     ,zpa,xpb,ypb,zpb,index,cg,tg,list,hlist,rlist,nb,work,iret
     &     ,errmsg)
      IF(iret .EQ. 1) RETURN

      IF(slv_redef) THEN
         CALL redefine_slv(ss_point,m1+1,ss_index,prsymb,mend,nbun
     &        ,slv_redef_type,mres,grppt,iret,errmsg)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
