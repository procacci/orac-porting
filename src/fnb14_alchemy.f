      SUBROUTINE fnb14_alchemy(ss_index,xp0,yp0,zp0,charge,nato,ecc12
     &     ,ecc6,cut,ewald,alphal,int14,int14p,type14,fudge,lj_fudge,fpx
     &     ,fpy,fpz,uconf_slt,uconf_slv,ucoul_slt,ucoul_slv,lambda
     &     ,lambdaq,lambda0,lambdaq0,dwrk,itask)

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
! Alchemic version.  Interactions within alchemical molecule 
! are normal (full) while the molecule is switched on/off. lambda scaled
! interactions occur only when one of the two interacting atoms is 
! alchemical. N.B. for the electrostatic part, the recpr 1-4 interactions
! as well as all the intra alchemical recpr interactions done with
! growing or decreasing charges need to 
! be subtraced out using a ferrf like scheme. This is done in mts_furier  
! where the ferrf is extended also to all intra-molecule  ij
! interactions with |i-j|>=4 

*======================= DECLARATIONS ==================================

      use unit
#ifdef _OMP_
      use omp_integr, only:array2_omp,nthr2
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,int14p,int14(2,*),type14(*),ss_index(*),itask
      REAL*8  xp0(nato),yp0(nato),zp0(nato),fpx(nato),fpy(nato),
     x        fpz(nato),ecc6(*),ecc12(*),alphal,charge(nato),
     x     fudge,lj_fudge,ucoul_slt,ucoul_slv,uconf_slt,uconf_slv
     &     ,lambda(*),lambdaq(*),lambda0(*),lambdaq0(*),dwrk
      LOGICAL ewald,cut

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,i2,j,lij,type
      REAL*8 xpi,ypi,zpi,xa,ya,za,rsq,rsp,rsqi,qforce
      REAL*8 xpj,ypj,zpj,ssvir,r6,r12,chrgei,chrgej
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,elj,ucon,ucou,ucoul(2),uconf(2)
      REAL*8 ucoul0(2),uconf0(2),furpar0,alch_fact,alch_fact0,alch_factq
     &     ,alch_factq0,ucon0,ucou0,aux,alch_factq2,ucoul20(2),furpar2
     &     ,furpar20,aux2,erff,alch_factq20,ucou2,ucou20
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         uconf(i)=0.0D+0
         ucoul(i)=0.0D+0
         uconf0(i)=0.0D+0
         ucoul0(i)=0.0D+0
      END DO
#ifndef _OMP_
      dwrk=0.d0     ! zeroing is outside parallel region in calling routine
#endif
      uconf_slt=0.0D0
      uconf_slv=0.0D0
      ucoul_slt=0.0D0
      ucoul_slv=0.0D0

      IF(int14p.EQ.0) RETURN
      twrtpi=2.0d0/DSQRT(pi)

      IF(.NOT.ewald.OR.(.NOT.cut)) THEN

*=======================================================================
*----- Do not take Ewald -----------------------------------------------
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
!           LENNARD_JONES
            IF(lambda(i1).eq.0.d0.and.lambda(i2).ne.0.d0) THEN 
              alch_fact=1.d0-dabs(lambda(i2))
              alch_fact0=1.d0-dabs(lambda0(i2))
            ELSE if(lambda(i1).ne.0.d0.and.lambda(i2).eq.0.d0) THEN  
              alch_fact=1.d0-dabs(lambda(i1))
              alch_fact0=1.d0-dabs(lambda0(i1))
!             normal potential if lambda's are zero 
            ELSE if (lambda(i1).eq.0.d0.and.lambda(i2).eq.0.d0) THEN   
              alch_fact=1.d0
              alch_fact0=1.d0
            ELSE if ((lambda(i1).ne.0.d0.and.lambda(i2).ne.0.d0).  
     &             AND.(lambda(i1).ne.lambda(i2))) THEN     
              alch_fact=0.d0 
              alch_fact0=0.d0 
!             on atoms (or off atoms) interact via the normal potential
            ELSE if ((lambda(i1).ne.0.d0.and.lambda(i2).ne.0.d0).  
     &             AND.(lambda(i1).eq.lambda(i2))) THEN   
              alch_fact=1.d0
              alch_fact0=1.d0
            END IF

!           CHARGES
            if(lambdaq(i1).eq.0.d0.and.lambdaq(i2).ne.0.d0) THEN 
              alch_factq=1.d0-dabs(lambdaq(i2))
              alch_factq0=1.d0-dabs(lambdaq0(i2))
            ELSE if(lambdaq(i1).ne.0.d0.and.lambdaq(i2).eq.0.d0) THEN  
              alch_factq=1.d0-dabs(lambdaq(i1))
              alch_factq0=1.d0-dabs(lambdaq0(i1))
!             normal potential if lambda's are zero 
            ELSE if (lambdaq(i1).eq.0.d0.and.lambdaq(i2).eq.0.d0) THEN   
              alch_factq=1.d0
              alch_factq0=1.d0
            ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(i2).ne.0.d0).  
     &             AND.(lambdaq(i1).ne.lambdaq(i2))) THEN     
              alch_factq=0.d0 
              alch_factq0=0.d0 
!             on atoms (or off atoms) interact via the normal potential
            ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(i2).ne.0.d0).  
     &             AND.(lambdaq(i1).eq.lambdaq(i2))) THEN   
              alch_factq=1.d0
              alch_factq0=1.d0
            END IF


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
            uconf(type)=uconf(type)+alch_fact*elj
            uconf0(type)=uconf0(type)+alch_fact*elj
            ssvir=fudge*chrgei*chrgej/rsp
            ucoul(type)=ucoul(type)+alch_factq*ssvir
            ucoul0(type)=ucoul0(type)+alch_factq0*ssvir
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
      ELSE IF(ewald) THEN

*=======================================================================
*----- Take Ewald for the electrostatic interaction --------------------
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
!           alch_fact controls the normal erfc part between qQ or qS
!           alch_fact2 sis used to subtract the (QQ-qq)*erf/r term 
            
!           LENNARD_JONES
            IF(lambda(i1).eq.0.d0.and.lambda(i2).ne.0.d0) THEN 
              alch_fact=1.d0-dabs(lambda(i2))
              alch_fact0=1.d0-dabs(lambda0(i2))
            ELSE if(lambda(i1).ne.0.d0.and.lambda(i2).eq.0.d0) THEN  
              alch_fact=1.d0-dabs(lambda(i1))
              alch_fact0=1.d0-dabs(lambda0(i1))
!             normal potential if lambda's are zero 
            ELSE if (lambda(i1).eq.0.d0.and.lambda(i2).eq.0.d0) THEN   
              alch_fact=1.d0
              alch_fact0=1.d0
            ELSE if ((lambda(i1).ne.0.d0.and.lambda(i2).ne.0.d0).  
     &             AND.(lambda(i1).ne.lambda(i2))) THEN     
              alch_fact=0.d0
              alch_fact0=0.d0
            ELSE if ((lambda(i1).ne.0.d0.and.lambda(i2).ne.0.d0).  
     &             AND.(lambda(i1).eq.lambda(i2))) THEN   
              alch_fact=1.d0
              alch_fact0=1.d0
            END IF

!           CHARGES
            if(lambdaq(i1).eq.0.d0.and.lambdaq(i2).ne.0.d0) THEN 
              alch_factq=1.d0-dabs(lambdaq(i2))
              alch_factq0=1.d0-dabs(lambdaq0(i2))
              alch_factq2=0.d0
              alch_factq20=0.0
            ELSE if(lambdaq(i1).ne.0.d0.and.lambdaq(i2).eq.0.d0) THEN  
              alch_factq=1.d0-dabs(lambdaq(i1))
              alch_factq0=1.d0-dabs(lambdaq0(i1))
              alch_factq2=0.d0
              alch_factq20=0.0
!             normal potential if lambda's are zero 
            ELSE if (lambdaq(i1).eq.0.d0.and.lambdaq(i2).eq.0.d0) THEN   
              alch_factq=1.d0
              alch_factq0=1.d0
              alch_factq2=0.d0
              alch_factq20=0.0
            ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(i2).ne.0.d0).  
     &             AND.(lambdaq(i1).ne.lambdaq(i2))) THEN     
              alch_factq=0.d0 
              alch_factq0=0.d0 
              alch_factq2=- (1.d0-dabs(lambdaq(i1)))*(1.
     &             -dabs(lambdaq(i2)))
              alch_factq20=- (1.d0-dabs(lambdaq0(i1)))*(1.
     &             -dabs(lambdaq0(i2)))
!           
!             on atoms (or off atoms) interact via the normal potential
            ELSE if ((lambdaq(i1).ne.0.d0.and.lambdaq(i2).ne.0.d0).  
     &             AND.(lambdaq(i1).eq.lambdaq(i2))) THEN   
              alch_factq=1.d0
              alch_factq0=1.d0
              alch_factq2=1.d0 - (1.d0-dabs(lambdaq(i1)))*(1.
     &             -dabs(lambdaq(i2)))
              alch_factq20=1.d0 - (1.d0-dabs(lambdaq0(i1)))*(1.
     &             -dabs(lambdaq0(i2)))
            END IF

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
            ssvir=alch_fact*12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
            qforce=lj_fudge*ssvir*rsqi
            aux= lj_fudge*(ecc12(lij)*r12-ecc6(lij)*r6)
            ucon=alch_fact*aux
            uconf(type)=uconf(type)+ucon
            ucon0=alch_fact0*aux
            uconf0(type)=uconf0(type)+ucon0
            alphar=alphal*rsp
            qt=1.0d0/(1.0d0+qp*alphar)
            expcst=dexp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
            furpar=alch_factq*fudge*chrgei*chrgej
            ucou=furpar*erfcst/rsp
            furpar0=alch_factq0*fudge*chrgei*chrgej
            ucou0=furpar0*erfcst/rsp
            ucoul(type)=ucoul(type)+ucou
            ucoul0(type)=ucoul0(type)+ucou0
            erff=1.d0-erfcst
            furpar2=alch_factq2*fudge*chrgei*chrgej
            ucou2=furpar2*erff/rsp
            ucoul(type)=ucoul(type)+ucou2
            furpar20=alch_factq20*fudge*chrgei*chrgej
            ucou20=furpar20*erff/rsp
            ucoul0(type)=ucoul0(type)+ucou20
            aux=(erfcst+twrtpi*alphar*expcst)*rspqi
            aux2=(erff-twrtpi*alphar*expcst)*rspqi
            qforce=qforce+furpar*aux + furpar2*aux2
            fpx(i1)=fpx(i1)+qforce*xa
            fpy(i1)=fpy(i1)+qforce*ya
            fpz(i1)=fpz(i1)+qforce*za
            fpx(i2)=fpx(i2)-qforce*xa
            fpy(i2)=fpy(i2)-qforce*ya
            fpz(i2)=fpz(i2)-qforce*za
20       CONTINUE
#ifdef _OMP_
          ENDIF
#endif
        END DO
      END IF

      uconf_slt=uconf(1)
      uconf_slv=uconf(2)
      ucoul_slt=ucoul(1)
      ucoul_slv=ucoul(2)
      do i=1,2
        dwrk=dwrk+uconf(i)+ucoul(i) - (uconf0(i)+ucoul0(i))
      end do

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
