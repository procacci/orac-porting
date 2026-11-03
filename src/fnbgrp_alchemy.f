      SUBROUTINE fnbgrp_alchemy(ss_index,xp0,yp0,zp0,charge,type,ecc12
     &     ,ecc6,ewald,alphal,ingrp,ingrpp,fpx,fpy,fpz,uconf_slt
     &     ,uconf_slv,ucoul_slt,ucoul_slv,lambdaq,itask)

************************************************************************
*                                                                      *
*     Compute the contribution from third neighbour (1-4) non          *
*     bonded interactions to forces and energies.                      *
*                                                                      *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     CHARGE  :  List of atomic charges for the solute.           (I)  *
*                >> real*8 CHARGE(NATO) <<                             *
*     CO      :  Transformation matrix from box coordinates       (I)  *
*                to orthogonal frame.                                  *
*                >> real*8 CO(3,3) <<                                  *
*     ECC12   :  List of L-J repulsive parameters.                (I)  *
*                >> real*8 ECC12(*) <<                                 *
*     ECC6    :  List of L-J attractive parameters.               (I)  *
*                >> real*8 ECC6(*) <<                                  *
*     CUT     :  Logical parameter. If .FALSE. all the non bonded (I)  *
*                interactions are included.                            *
*                >> logical*4 CUT <<                                   *
*     EWALD   :  Logical parameter. If .TRUE. the electrostatic   (I)  *
*                interaction is compute with Ewald.                    *
*                >> logical*4 EWALD <<                                 *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     FPX     :  Forces for each atom of the macromolecule.      (I/O) *
*     FPY        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*     FPZ                                                              *
*                                                                      *
*     UCONF   :  Configurational energy.                           (O) *
*     UCOUL   :  Coulombic energy.                                 (O) *
*                                                                      *
*---- Last update 05/22/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
#ifdef _OMP_
      use omp_integr, only:nthr2,array2_omp
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,ingrpp,type(*),ingrp(2,*),ss_index(*),itask
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),ecc6(*),ecc12(*)
     &     ,lambdaq(*),alphal,charge(*),uconf_slt,ucoul_slt,uconf_slv
     &     ,ucoul_slv
      LOGICAL ewald

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,i2,j,lij,li,lj,typei
      REAL*8 xpi,ypi,zpi,xa,ya,za,rsq,rsp,rsqi,qforce,alphal0
      REAL*8 xpj,ypj,zpj,ssvir,r6,r12,chrgei,chrgej
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,ucoul(2),uconf(2)
      REAL*8 rspqi,alphar,furpar,twrtpi,elj,ucon,ucou
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
        uconf(i)=0.0D+0
        ucoul(i)=0.0D+0
      END DO
      uconf_slt=0.0D0
      uconf_slv=0.0D0
      ucoul_slt=0.0D0
      ucoul_slv=0.0D0
      IF(ingrpp.EQ.0) RETURN
      twrtpi=2.0d0/DSQRT(pi)

      IF(.NOT.ewald) THEN

*=======================================================================
*----- Do not take Ewald -----------------------------------------------
*=======================================================================

        DO  i=1,ingrpp
#ifdef _OMP_
          if(array2_omp(i,itask)) THEN 
#endif
            i1=ingrp(1,i)
            i2=ingrp(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            typei=ss_index(i1)
            lj=MAX0(type(i1),type(i2))
            li=MIN0(type(i1),type(i2))
            lij=lj*(lj-1)/2+li
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            xpj=xp0(i2)
            ypj=yp0(i2)
            zpj=zp0(i2)
            chrgei=charge(i1)
            chrgej=charge(i2)
            xa=xpi-xpj
            ya=ypi-ypj
            za=zpi-zpj
            rsq=xa**2+ya**2+za**2
            rsp=DSQRT(rsq)
            rsqi=1.0d0/rsq
            r6=rsqi*rsqi*rsqi
            r12=r6*r6
            ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
            qforce=ssvir*rsqi
            elj=ecc12(lij)*r12-ecc6(lij)*r6
            uconf(typei)=uconf(typei)+elj
            ssvir=chrgei*chrgej/rsp
            ucoul(typei)=ucoul(typei)+ssvir
            qforce=qforce+ssvir*rsqi
            fpx(i1)=fpx(i1)+qforce*xa
            fpy(i1)=fpy(i1)+qforce*ya
            fpz(i1)=fpz(i1)+qforce*za
            fpx(i2)=fpx(i2)-qforce*xa
            fpy(i2)=fpy(i2)-qforce*ya
            fpz(i2)=fpz(i2)-qforce*za
10          CONTINUE
#ifdef _OMP_
          ENDIF
#endif        
        END DO
      ELSE

*=======================================================================
*----- Take Ewald for the electrostatic interaction --------------------
*=======================================================================

        DO  i=1,ingrpp
#ifdef _OMP_
          if(array2_omp(i,itask)) THEN 
#endif
            i1=ingrp(1,i)
            i2=ingrp(2,i)
            
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            if(dabs(lambdaq(i1)).eq.1.d0.and.
     &           dabs(lambdaq(i2)).eq.1.d0) THEN 
              alphal0=0.d0
            ELSE
              alphal0=alphal
            END IF
            typei=ss_index(i1)
            lj=MAX0(type(i1),type(i2))
            li=MIN0(type(i1),type(i2))
            lij=lj*(lj-1)/2+li
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            xpj=xp0(i2)
            ypj=yp0(i2)
            zpj=zp0(i2)
            chrgei=charge(i1)
            chrgej=charge(i2)
            xa=xpi-xpj
            ya=ypi-ypj
            za=zpi-zpj
            rsq=xa**2+ya**2+za**2
            rsp=DSQRT(rsq)
            rsqi=1.0d0/rsq
            rspqi=rsqi/rsp
            r6=rsqi*rsqi*rsqi
            r12=r6*r6
            ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
            qforce=ssvir*rsqi
            ucon=ecc12(lij)*r12-ecc6(lij)*r6
            uconf(typei)=uconf(typei)+ucon
            alphar=alphal*rsp
            qt=1.0d0/(1.0d0+qp*alphar)
            expcst=dexp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x           *qt+a1)*qt*expcst
            furpar=chrgei*chrgej
            ucou=furpar*erfcst/rsp
            ucoul(typei)=ucoul(typei)+ucou
            qforce=qforce+furpar*(erfcst+twrtpi*alphar
     x           *expcst)*rspqi
            fpx(i1)=fpx(i1)+qforce*xa
            fpy(i1)=fpy(i1)+qforce*ya
            fpz(i1)=fpz(i1)+qforce*za
            fpx(i2)=fpx(i2)-qforce*xa
            fpy(i2)=fpy(i2)-qforce*ya
            fpz(i2)=fpz(i2)-qforce*za
20          CONTINUE
#ifdef _OMP_
          ENDIF
#endif
        END DO
      END IF
      
      uconf_slt=uconf(1)
      uconf_slv=uconf(2)
      ucoul_slt=ucoul(1)
      ucoul_slv=ucoul(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
