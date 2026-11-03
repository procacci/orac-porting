      SUBROUTINE fnb14_rem(ss_index,xp0,yp0,zp0,charge,nato,ecc12,ecc6
     &     ,cut,ewald,alphal,int14,int14p,type14,fudge,lj_fudge,fpx,fpy
     &     ,fpz,uconf_slt,uconf_slv,ucoul_slt,ucoul_slv,remint14,itask
     &     ,seg_ene)

************************************************************************
*                                                                      *
*     Compute the contribution from third neighbour (1-4) non          *
*     bonded interactions to forces and energies. (REM version)        *
*                                                                      *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     CHARGE  :  List of atomic charges for the solute.           (I)  *
*                >> real*8 CHARGE(NATO) <<                             *
*     NATO    :  Number of atoms forming the solute.              (I)  *
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
*     INT14   :  List of all 1-4 interactions.                    (I)  *
*                >> integer INT14(2,*) <<                              *
*     INT14P  :  Number of 1-4 interactions.                      (I)  *
*     TYPE14  :  List of interaction type for each 1-4 interaction.(I) *
*                >> integer*4 TYPE14(*) <<                             *
*     FUDGE   :  Fudge parameter. It multiplies the Lennard-Jones (I)  *
*                interaction.                                          *
*                >> real*8 FUDGE <<                                    *
*     REMINT14:  Logical array. If .TRUE. the 14 interaction for the   *
*                pair is multiplied by the appropriate REM             *
*                scaling factor                                        *
*                >> LOGICAL REMINT(M3) <<                              *
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
      USE rem
#ifdef _OMP_
      use omp_integr, only:nthr2,array2_omp
#endif
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,int14p,int14(2,*),type14(*),ss_index(*),itask
      REAL*8  xp0(nato),yp0(nato),zp0(nato),fpx(nato),fpy(nato),
     x        fpz(nato),ecc6(*),ecc12(*),alphal,charge(nato),
     x     fudge,lj_fudge,ucoul_slt,ucoul_slv,uconf_slt
     &     ,uconf_slv,seg_ene
      LOGICAL ewald,cut
      LOGICAL remint14(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,i2,j,lij,type
      REAL*8 xpi,ypi,zpi,xa,ya,za,rsq,rsp,rsqi,qforce
      REAL*8 xpj,ypj,zpj,ssvir,r6,r12,chrgei,chrgej
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,elj,ucon,ucou,ucoul(2),uconf(2)
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

      IF(int14p.EQ.0) RETURN
      twrtpi=2.0d0/DSQRT(pi)

      seg_ene=0.d0
      IF(.NOT.ewald.OR.(.NOT.cut)) THEN

*=======================================================================
*----- Do not take Ewald -----------------------------------------------
*=======================================================================

         DO i=1,int14p
#ifdef _OMP_
           if(array2_omp(i,itask)) THEN 
#endif
            i1=int14(1,i)
            i2=int14(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            type=ss_index(i1)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            xpj=xp0(i2)
            ypj=yp0(i2)
            zpj=zp0(i2)
            lij=type14(i)
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
            qforce=lj_fudge*ssvir*rsqi
            elj=lj_fudge*(ecc12(lij)*r12-ecc6(lij)*r6)
            uconf(type)=uconf(type)+elj
            ssvir=fudge*chrgei*chrgej/rsp
            ucoul(type)=ucoul(type)+ssvir
            qforce=qforce+ssvir*rsqi
            if(remint14(i)) then
              seg_ene = seg_ene + elj + ssvir
              qforce = qforce * rem_factor(2)
            endif
            fpx(i1)=fpx(i1)+qforce*xa
            fpy(i1)=fpy(i1)+qforce*ya
            fpz(i1)=fpz(i1)+qforce*za
            fpx(i2)=fpx(i2)-qforce*xa
            fpy(i2)=fpy(i2)-qforce*ya
            fpz(i2)=fpz(i2)-qforce*za
#ifdef _OMP_
          ENDIF
#endif
10      CONTINUE
       END DO
      ELSE IF(ewald) THEN

*=======================================================================
*----- Take Ewald for the electrostatic interaction --------------------
*=======================================================================

         DO  i=1,int14p
#ifdef _OMP_
           if(array2_omp(i,itask)) THEN 
#endif
            i1=int14(1,i)
            i2=int14(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            type=ss_index(i1)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            xpj=xp0(i2)
            ypj=yp0(i2)
            zpj=zp0(i2)
            lij=type14(i)
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
            qforce=lj_fudge*ssvir*rsqi
            ucon=lj_fudge*(ecc12(lij)*r12-ecc6(lij)*r6)
            uconf(type)=uconf(type)+ucon
            alphar=alphal*rsp
            qt=1.0d0/(1.0d0+qp*alphar)
            expcst=dexp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
            furpar=fudge*chrgei*chrgej
            ucou=furpar*erfcst/rsp
            ucoul(type)=ucoul(type)+ucou
            qforce=qforce+furpar*(erfcst+twrtpi*alphar*expcst)*rspqi
            if(remint14(i)) then
              seg_ene = seg_ene + ucon + ucou
              qforce = qforce * rem_factor(2)
            endif
            fpx(i1)=fpx(i1)+qforce*xa
            fpy(i1)=fpy(i1)+qforce*ya
            fpz(i1)=fpz(i1)+qforce*za
            fpx(i2)=fpx(i2)-qforce*xa
            fpy(i2)=fpy(i2)-qforce*ya
            fpz(i2)=fpz(i2)-qforce*za
#ifdef _OMP_
          ENDIF
#endif
20        CONTINUE
        END DO
      END IF

      uconf_slt=uconf(1)
      uconf_slv=uconf(2)
      ucoul_slt=ucoul(1)
      ucoul_slv=ucoul(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
