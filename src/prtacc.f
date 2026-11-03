      SUBROUTINE prtacc(pucek,puhyd,puconf,pueng,pucoul,sftalp,fsbond
     &     ,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond,uitors,uptors
     &     ,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal,ucek
     &     ,upconf,upcoul,uslvbon,uslvben,uslvtor,uslvitor,ugyr,rms,temp
     &     ,temph,tcm,rcm,tempt,tempr,temppr,gra,ucepr,prt
     &     ,press,pressc,pressk,temppra,errca,errhe,errbc,erral,drpca
     &     ,drpbc,drphe,drpal,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot
     &     ,sum_eslvint,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors
     &     ,sum_tote,sum_ucek,sum_temp,sum_tempt,sum_tempr,sum_temppr
     &     ,sum_temph,sum_rms,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul
     &     ,sum_percip,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor
     &     ,sum_pnbd,sum_pebnd,sum_pepot,sum_ptote,sum_gr,sum_epcoul
     &     ,sum_epconf,sum_co,sum_st,sum_presst,sum_press,sum_pressc
     &     ,sum_pressk,sum_volume,sum_pv,sum_temppra,ssm_econf,ssm_ecoul
     &     ,ssm_enbnd,ssm_etotpot,ssm_eslvint,ssm_eebond,ssm_eebend
     &     ,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek,ssm_temp,ssm_tempt
     &     ,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms,ssm_pecek,ssm_pehyd
     &     ,ssm_peconf,ssm_pecoul,ssm_percip,ssm_enb14,ssm_ebend
     &     ,ssm_ebond,ssm_eitor,ssm_eptor,ssm_pnbd,ssm_pebnd,ssm_pepot
     &     ,ssm_ptote,ssm_gr,ssm_epcoul,ssm_epconf,ssm_co,ssm_st
     &     ,ssm_presst,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &     ,ssm_temppra,energy,enisolv,mstep,nstep)


*----------------------- ARGUMENTS -------------------------------------

      use unit
      use parst
      use cpropar

      IMPLICIT none

      REAL*8 pucek,puhyd,puconf,pueng,pucoul,sftalp,fsbond,purecp,fsbend
     x     ,fsin14,unb14,cnb14,ubend,ubond,uitors,uptors,pubnd,uceh,hpot
     x     ,ucoul,uconf,urecp,ureal,fsrtal,ucek,upconf,upcoul,rms,temp
     x     ,temph,tcm,rcm,tempt,tempr,temppr,gra,ucepr,prt(3
     &     ,3),press,pressc,pressk,temppra,energy,uslvbon,uslvben
     &     ,uslvtor,uslvitor,ugyr,enisolv
      REAL*8 errca(*),errhe(*),errbc(*),erral(*),drpca(*),drpbc(*)
     &     ,drphe(*),drpal(*)

      REAL*8 sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_tote,sum_ucek
     x     ,sum_temp,sum_tempt,sum_tempr,sum_temppr,sum_temph,sum_rms
     x     ,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul,sum_percip
     x     ,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor,sum_pnbd
     x     ,sum_pebnd,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     x     ,sum_co(3,3),sum_st(3,3),sum_presst,sum_press,sum_pressc
     x     ,sum_pressk,sum_volume,sum_pv,sum_temppra,sum_eslvint
     &     ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors

      REAL*8 ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_tote,ssm_ucek
     x     ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
     x     ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     x     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor,ssm_pnbd
     x     ,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul,ssm_epconf
     x     ,ssm_co(3,3),ssm_st(3,3),ssm_presst,ssm_press,ssm_pressc
     x     ,ssm_pressk,ssm_volume,ssm_pv,ssm_temppra,ssm_eslvint
     &     ,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors
      INTEGER mstep,nstep,i,j

*-------------------- Local Variables  ---------------------------------

      INTEGER ninner,nprop1
      REAL*8  pecek,pehyd,pecoul,percip,enb14,ebend,ebond,eitor,eptor
     &     ,pnbd,pebnd,pepot,ptote,econf,ecoul,erecp,ereal,econd,tote
     &     ,epconf,epcoul,peconf,ssm,step,avg,FLUCT,fqrtal,fnstep,ecek
     &     ,pv,presst,eslvintra,egyr,mrject,eebond,eebend,eeptors
     &     ,eeitors,totpot,interm
      REAL*8 avg_econf,avg_ecoul,avg_enbnd,avg_etotpot,avg_tote,avg_ucek
     &     ,avg_temp,avg_tempt,avg_tempr,avg_temppr,avg_temph,avg_rmsp
     &     ,avg_pecek,avg_pehyd,avg_peconf,avg_pecoul,avg_percip
     &     ,avg_enb14,avg_ebend,avg_ebond,avg_eitor,avg_eptor,avg_pnbd
     &     ,avg_pebnd,avg_pepot,avg_ptote,avg_epcoul,avg_epconf,avg_co(3
     &     ,3),avg_st(3,3),avg_presst,avg_press,avg_pressc,avg_pressk
     &     ,avg_volume,avg_pv,avg_temppra,avg_eslvint,avg_eebond
     &     ,avg_eebend,avg_eeptors,avg_eeitors
      REAL*8 flc_econf,flc_ecoul,flc_enbnd,flc_etotpot,flc_tote,flc_ucek
     &     ,flc_temp,flc_tempt,flc_tempr,flc_temppr,flc_temph,flc_rms
     &     ,flc_pecek,flc_pehyd,flc_peconf,flc_pecoul,flc_percip
     &     ,flc_enb14,flc_ebend,flc_ebond,flc_eitor,flc_eptor,flc_pnbd
     &     ,flc_pebnd,flc_pepot,flc_ptote,flc_epcoul,flc_epconf,flc_co(3
     &     ,3),flc_st(3,3),flc_presst,flc_press,flc_pressc,flc_pressk
     &     ,flc_volume,flc_pv,flc_temppra,flc_eslvint,flc_eebond
     &     ,flc_eebend,flc_eeptors,flc_eeitors
      REAL*8  a,b,c,alpha,bbeta,gamma,fstep,gr,io,co2(3,3)
      DATA ninner/0/
      FLUCT(ssm,step,avg)=DSQRT(DABS((ssm-step*avg*avg)/step))

      IF(abmd_tors) THEN
         gr=DMOD(gra*180.0D0/pi,360.0D0)
      ELSE
         gr=gra
      END IF
      IF(md_respa) THEN
         mrject=nrject*lrespa*mrespa
         fstep=time*DFLOAT(mstep)/DFLOAT(mrespa*lrespa)
         nprop1=nprop*lrespa*mrespa
      ELSE
         mrject=nrject
         fstep=time*DFLOAT(mstep)
         nprop1=nprop
      END IF

      pv=0.0D0
      IF(pressure .OR. cpress .OR. lberendsen) THEN
         pv=(pext*volume*efact)/1000.0D0
         DO i=1,3
            DO j=1,3
               prt(i,j)=prt(i,j)*unitp/(3.0D0*volume)
            END DO
         END DO
      END IF
      
      IF(clewld) percip=(purecp+fsbond+sftalp+fsbend+fsin14)*efact
     &     /1000.0D0

      IF(MOD(mstep,nprint).EQ.0)THEN
         WRITE(kprint,'(/)')
      END IF
      IF(slt_exist.AND.(.NOT.slv_exist)) THEN

*=======================================================================
*==== WRITE INTERMEDIATE RESULTS WHEN REQUIRED =========================
*=======================================================================

         pecek=pucek/1000.0d0
         pehyd=puhyd*efact/1000.0d0
         peconf=puconf*efact/1000.0d0
         egyr=ugyr*efact/1000.0D0
         IF(uenerg) pueng=pueng*efact/1000.0D0
         pecoul=pucoul*efact/1000.0D0
         enb14=(unb14+cnb14)*efact/1000.0d0
         ebend=ubend*efact/1000.0D0
         ebond=ubond*efact/1000.0D0
         eitor=uitors*efact/1000.0D0
         eptor=uptors*efact/1000.0D0
         pnbd=peconf+pecoul+pehyd+pueng
         pebnd=pubnd*efact/1000.0d0
         pepot=pnbd+pebnd
         ptote=pnbd+pecek+pebnd
         IF(hoover .OR. thermos) THEN
            ptote=ptote+(uceh+hpot*efact)/1000.0d+0
         END IF
         IF(cpress .OR. lberendsen) THEN
            ptote=ptote+pv+ucepr/1000.0D0
         END IF
         if(agbnp.or.lsasa) THEN 
           ptote=ptote+enisolv*efact/1000.0D0
         ENDIF 
         IF(MOD(mstep,nprint).EQ.0)THEN
            IF(.NOT.abmd) THEN
               IF(.NOT.clewld) THEN
                  IF(hoover .OR. thermos) THEN
                     WRITE(kprint,12050) fstep,ptote,pepot,
     x                    pecoul,pnbd,enb14,pebnd,ebond,ebend,eitor,
     x                    eptor,temp,temph,rcm,tcm
                  ELSE
                    IF(.not.(agbnp.or.lsasa)) THEN
                      WRITE(kprint,12001) fstep,ptote
     &                     ,pepot,pecoul,pnbd,enb14,pebnd,ebond,ebend
     &                     ,eitor,eptor,temp,rcm,tcm
#ifndef _PGI_
                      IF(isnan(ptote)) THEN
                        write(6,*) ptote,
     &            "ptote is NaN. Check coordinates; Program stop."
                        STOP
                      ENDIF
#endif
                    ELSE
                      WRITE(kprint,12071) fstep,ptote,pepot,
     x                     pecoul,pnbd,enb14,pebnd,ebond,ebend,eitor,
     x                     eptor,temp,rcm,tcm,enisolv*efact/1000.d0
                    ENDIF
                  END IF
               ELSE
                  IF(hoover .OR. thermos) THEN
                     WRITE(kprint,12051) fstep,ptote,pepot,
     x                    pecoul,percip,pnbd,enb14,pebnd,ebond,ebend
     &                    ,eitor,eptor,temp,temph,rcm,tcm
                  ELSE
                     WRITE(kprint,12006) fstep,ptote,pepot,
     x                    pecoul,percip,pnbd,enb14,pebnd,ebond,ebend
     &                    ,eitor,eptor,temp,rcm,tcm
                  END IF
               END IF
            ELSE
               IF(hoover .OR. thermos) THEN
                  WRITE(kprint,12054) fstep,ptote,pepot,pecoul,pnbd
     x                 ,enb14,pebnd,ebond,ebend,eitor,eptor,temp,temph
     x                 ,rcm,tcm,gr,egyr
               ELSE
                  WRITE(kprint,12004) fstep,ptote,pepot,pecoul,pnbd
     x                 ,enb14,pebnd,ebond,ebend,eitor,eptor,temp,rcm
     x                 ,tcm,gr,egyr
               END IF
               IF(fold .AND. nabmd .NE. 0) THEN
                  WRITE(kabmd,80300) fstep,gr,egyr,pepot
               END IF
            END IF
         END IF
      END IF
      IF(anxrms) THEN
         ninner=ninner+1
         IF(MOD(mstep,nxrms).EQ. 0) THEN
            IF(anprot) THEN
               WRITE(kxrms,50000) fstep
               IF(anxca) CALL write_xrms(kxrms,annpro,'CA',errca)
               IF(anxbc) CALL write_xrms(kxrms,annpro,'BC',errbc)
               IF(anxhe) CALL write_xrms(kxrms,annpro,'HE',errhe)
               IF(anxal) CALL write_xrms(kxrms,annpro,'AL',erral)
            ELSE
               WRITE(kxrms,50000) fstep
               IF(anxca) CALL write_xrms(kxrms,nprot,'CA',errca)
               IF(anxbc) CALL write_xrms(kxrms,nprot,'BC',errbc)
               IF(anxhe) CALL write_xrms(kxrms,nprot,'HE',errhe)
               IF(anxal) CALL write_xrms(kxrms,nprot,'AL',erral)
            END IF
            REWIND kxrms_atm
            IF(anxca) CALL write_xrms_atm(kxrms_atm,ntap,'CA',drpca
     &           ,ninner,fstep,ngrp,grppt,nres(1,1))
            IF(anxbc) CALL write_xrms_atm(kxrms_atm,ntap,'BC',drpbc
     &           ,ninner,fstep,ngrp,grppt,nres(1,1))
            IF(anxhe) CALL write_xrms_atm(kxrms_atm,ntap,'HE',drphe
     &           ,ninner,fstep,ngrp,grppt,nres(1,1))
            IF(anxal) CALL write_xrms_atm(kxrms_atm,ntap,'AL',drpal
     &           ,ninner,fstep,ngrp,grppt,nres(1,1))
         END IF
      END IF
         
      IF(menerg .GT. 0) THEN
         IF(MOD(mstep,menerg).EQ.0.AND.mstep.GT.mrject)THEN
            WRITE(kenerg,'(7e15.7)') pepot,pnbd,pebnd,peconf,
     x              pecoul,pehyd,pueng
         END IF
      END IF
      IF(slv_exist.AND.(.NOT.slt_exist)) THEN
          econf=uconf*efact/1000.0D0
          ecoul=ucoul*efact/1000.0D0
          erecp=urecp*efact/1000.0D0
          ereal=ureal*efact/1000.0D0
          econd=econf-ecoul
          tote=econf+ucek/1000.0d+0
          ecek=ucek/1000.0D0
          egyr=ugyr*efact/1000.0D0
          IF(hoover .OR. thermos) THEN
              tote=tote+(uceh+hpot*efact)/1000.0D0
           END IF
           
          IF(cpress .OR. lberendsen) THEN
             tote=tote+pv+ucepr/1000.0D0
          END IF
          tote=tote+(uslvbon+uslvben+uslvtor+uslvitor)*efact/1000.0D0
          eslvintra=(uslvbon+uslvben+uslvtor+uslvitor)*efact/1000.0D0
          eebond=uslvbon*efact/1000.0D0
          eebend=uslvben*efact/1000.0D0
          eeptors=uslvtor*efact/1000.0D0
          eeitors=uslvitor*efact/1000.0D0
          totpot=eslvintra+econf
c$$$          WRITE(6,'(f30.25)') totpot,ecek
c$$$          WRITE(6,'(3f25.20)') ((prt(i,j)/unitp,j=1,3),i=1,3)

          IF(.NOT.abmd) THEN
             IF(MOD(mstep,nprint).EQ.0)THEN
                IF(.NOT. (hoover .OR. thermos)) THEN
                   WRITE(kprint,12000) fstep,tote,totpot,econf,ecoul
     &                  ,eslvintra,eebond,eebend,eeptors,eeitors,tempt
     &                  ,tempr,temp,rms
                ELSE
                   WRITE(kprint,12010) fstep,tote,totpot,econf,ecoul
     &                  ,eslvintra,eebond,eebend,eeptors,eeitors,tempt
     &                  ,tempr,temp,temph,rms
                END IF
             END IF
          ELSE
             IF(MOD(mstep,nprint).EQ.0)THEN
                IF(.NOT. (hoover .OR. thermos)) THEN
                   WRITE(kprint,13000) fstep,tote,totpot,econf,ecoul
     &                  ,eslvintra,eebond,eebend,eeptors,eeitors,tempt
     &                  ,tempr,temp,gr,egyr
                ELSE
                   WRITE(kprint,13010) fstep,tote,totpot,econf,ecoul
     &                  ,eslvintra,eebond,eebend,eeptors,eeitors,tempt
     &                  ,tempr,temp,temph,gr,egyr
                END IF
             END IF
          END IF
          IF(abmd_tors .AND. nabmd .NE. 0) THEN
             IF(MOD(mstep,nabmd) .EQ. 0) THEN
                WRITE(kabmd,80200) fstep,gra*180.0D0/pi,egyr,econf
     &               ,eeptors
             END IF
          END IF
      END IF

      IF(slv_exist.AND.slt_exist) THEN

*=======================================================================
*====== Adds up all the solute-solute contributions ====================
*=======================================================================
         
         pecek=pucek/1000.0d0
         pehyd=puhyd*efact/1000.0d0
         peconf=puconf*efact/1000.0d0
         pecoul=pucoul*efact/1000.0d0
         egyr=ugyr*efact/1000.0D0
         enb14=(unb14+cnb14)*efact/1000.0d0
         ebend=ubend*efact/1000.0D0
         ebond=ubond*efact/1000.0D0
         eitor=uitors*efact/1000.0D0
         eptor=uptors*efact/1000.0D0
         pnbd=peconf+pecoul+pehyd
         pebnd=pubnd*efact/1000.0d0
         pepot=pnbd+pebnd
         ptote=pnbd+pecek+pebnd
         if(agbnp.or.lsasa) THEN 
           ptote=ptote+enisolv*efact/1000.0D0
         ENDIF 

*=======================================================================
*====== Adds up all the slv_existt-slv_existt contributions ==================
*=======================================================================

         econf=uconf*efact/1000.0d+0
         ecoul=ucoul*efact/1000.0d+0
         erecp=urecp*efact/1000.0d+0
         ereal=ureal*efact/1000.0d+0

*=======================================================================
*====== Adds up all the slv_existt-solute contributions ===================
*=======================================================================

         epconf=upconf*efact/1000.0D0
         epcoul=upcoul*efact/1000.0D0

*=======================================================================
*====== Finally calculate the total energy =============================
*=======================================================================

         tote=econf+ucek/1000.0D0+ptote+epconf
         ecek=ucek/1000.0D0
         IF(hoover .OR. thermos) THEN
            tote=tote+(uceh+hpot*efact)/1000.0d+0
         END IF
         IF(cpress .OR. lberendsen) THEN
            tote=tote+(pv+ucepr/1000.0D0)
         END IF
         eslvintra=(uslvbon+uslvben+uslvtor+uslvitor)*efact/1000.0D0
         tote=tote+eslvintra
         interm=econf
         econf=econf+eslvintra
         eebond=uslvbon*efact/1000.0D0
         eebend=uslvben*efact/1000.0D0
         eeptors=uslvtor*efact/1000.0D0
         eeitors=uslvitor*efact/1000.0D0
c$$$         CALL dgemm('N','T',3,3,3,1.0D0,prt,3,oc,3,0.0D0
c$$$     &        ,co2,3)
c$$$         WRITE(6,'(3f20.10)') ((co2(i,j)/unitp,j=1,3),i=1,3)
c$$$         WRITE(6,'(f30.25)') tote*1000.0D0/efact,pebnd*1000.0D0/efact
c$$$     &        ,eslvintra*1000.0D0/efact
c$$$         WRITE(6,'(3f25.20)') ((prt(i,j)/unitp,j=1,3),i=1,3)
         IF(MOD(mstep,nprint).EQ.0)THEN
            IF(.NOT. (hoover .OR. thermos)) THEN
              IF(.not.(agbnp.or.lsasa)) THEN 
                WRITE(kprint,12002) fstep,tote,econf,interm,ecoul
     &               ,eslvintra,eebond,eebend,eeptors,eeitors,ptote
     &               ,pepot,pecoul,peconf,pehyd,pebnd,ebond,ebend,eitor
     &               ,eptor,epconf,epcoul,temppr,tempt,tempr,temp
              ELSE
                WRITE(kprint,12022) fstep,tote,econf,interm,ecoul
     &               ,eslvintra,eebond,eebend,eeptors,eeitors,ptote
     &               ,pepot,pecoul,peconf,pehyd,pebnd,ebond,ebend,eitor
     &               ,eptor,epconf,epcoul,temppr,tempt,tempr,temp
     &               ,enisolv*efact/1000.d0
              end if
            ELSE
               WRITE(kprint,12003) fstep,tote,econf,interm,ecoul
     &              ,eslvintra,eebond,eebend,eeptors,eeitors,ptote,pepot
     &              ,pecoul,peconf,pehyd,pebnd,ebond,ebend,eitor,eptor
     &              ,epconf,epcoul,temppr,tempt,tempr,temp,temph
            END IF
         END IF
      END IF
       
*=======================================================================
*==== Calculate stress tensor and pressure when required ===============
*=======================================================================

      IF(cpress .OR. lberendsen .OR. pressure) THEN
         CALL rotb(a,b,c,alpha,bbeta,gamma,co)
         presst=press+pressk
         IF(MOD(mstep,nprint).EQ.0)THEN

            WRITE(kprint,12200) presst,press,pressk,temppra,volume,pv
            WRITE(kprint,12300) a,b,c,prt(1,1),prt(1,2),prt(1,3),alpha
     &           ,bbeta,gamma,prt(2,1),prt(2,2),prt(2,3),prt(3,1),prt(3
     &           ,2),prt(3,3)
         END IF
      END IF

*=======================================================================
*==== ACCUMULATE AVERAGES AND WRITE SUBAVERAGES WHEN REQUIRED ==========
*=======================================================================

      IF(mstep-mrject .GT. 0 .AND. mdsim)THEN
         IF(slt_exist .AND. (.NOT. slv_exist) ) THEN
            sum_pecek=sum_pecek+pecek
            sum_pehyd=sum_pehyd+pehyd
            sum_peconf=sum_peconf+peconf
            sum_pecoul=sum_pecoul+pecoul
            sum_percip=sum_percip+percip
            sum_enb14=sum_enb14+enb14
            sum_ebend=sum_ebend+ebend
            sum_ebond=sum_ebond+ebond
            sum_eitor=sum_eitor+eitor
            sum_eptor=sum_eptor+eptor
            sum_pnbd=sum_pnbd+pnbd
            sum_pebnd=sum_pebnd+pebnd
            sum_pepot=sum_pepot+pepot
            sum_ptote=sum_ptote+ptote
            sum_temp=sum_temp+temp
            sum_temph=sum_temph+temph
            sum_gr=sum_gr+gr
            ssm_pecek=ssm_pecek+pecek**2
            ssm_pehyd=ssm_pehyd+pehyd**2
            ssm_peconf=ssm_peconf+peconf**2
            ssm_pecoul=ssm_pecoul+pecoul**2
            ssm_percip=ssm_percip+percip**2
            ssm_enb14=ssm_enb14+enb14**2
            ssm_ebend=ssm_ebend+ebend**2
            ssm_ebond=ssm_ebond+ebond**2
            ssm_eitor=ssm_eitor+eitor**2
            ssm_eptor=ssm_eptor+eptor**2
            ssm_pnbd=ssm_pnbd+pnbd**2
            ssm_pebnd=ssm_pebnd+pebnd**2
            ssm_pepot=ssm_pepot+pepot**2
            ssm_ptote=ssm_ptote+ptote**2
            ssm_temp=ssm_temp+temp**2
            ssm_temph=ssm_temph+temph**2
            ssm_gr=ssm_gr+gr**2
            
            IF(MOD(mstep,nprop1).EQ.0)THEN
               fnstep=DFLOAT(mstep)
               avg_temp=sum_temp/fnstep
               avg_temph=sum_temph/fnstep
               avg_pecek=sum_pecek/fnstep
               avg_pehyd=sum_pehyd/fnstep
               avg_peconf=sum_peconf/fnstep
               avg_pecoul=sum_pecoul/fnstep
               avg_percip=sum_percip/fnstep
               avg_enb14=sum_enb14/fnstep
               avg_ebend=sum_ebend/fnstep
               avg_ebond=sum_ebond/fnstep
               avg_eitor=sum_eitor/fnstep
               avg_eptor=sum_eptor/fnstep
               avg_pnbd=sum_pnbd/fnstep
               avg_pebnd=sum_pebnd/fnstep
               avg_pepot=sum_pepot/fnstep
               avg_ptote=sum_ptote/fnstep
               
               flc_temp=FLUCT(ssm_temp,fnstep,avg_temp)
               flc_temph=FLUCT(ssm_temph,fnstep,avg_temph)
               flc_pecek=FLUCT(ssm_pecek,fnstep ,avg_pecek)
               flc_pehyd=FLUCT(ssm_pehyd,fnstep ,avg_pehyd)
               flc_peconf=FLUCT(ssm_peconf,fnstep ,avg_peconf)
               flc_pecoul=FLUCT(ssm_pecoul,fnstep ,avg_pecoul)
               flc_percip=FLUCT(ssm_percip,fnstep ,avg_percip)
               flc_enb14=FLUCT(ssm_enb14,fnstep ,avg_enb14)
               flc_ebend=FLUCT(ssm_ebend,fnstep ,avg_ebend)
               flc_ebond=FLUCT(ssm_ebond,fnstep ,avg_ebond)
               flc_eitor=FLUCT(ssm_eitor,fnstep ,avg_eitor)
               flc_eptor=FLUCT(ssm_eptor,fnstep ,avg_eptor)
               flc_pnbd=FLUCT(ssm_pnbd,fnstep ,avg_pnbd)
               flc_pebnd=FLUCT(ssm_pebnd,fnstep ,avg_pebnd)
               flc_pepot=FLUCT(ssm_pepot,fnstep ,avg_pepot)
               flc_ptote=FLUCT(ssm_ptote,fnstep ,avg_ptote)
               WRITE(kprint,40000) fstep
               WRITE(kprint,40100) 
               IF(.NOT. (hoover  .OR. thermos)) THEN
                  WRITE(kprint,22001) avg_ptote,flc_ptote,
     x                 avg_pepot,flc_pepot,avg_pecoul,flc_pecoul,
     x                 avg_pnbd,flc_pnbd,avg_enb14,flc_enb14,
     x                 avg_pebnd,flc_pebnd,avg_ebond,flc_ebond,
     x                 avg_ebend,flc_ebend,avg_eitor,flc_eitor,
     x                 avg_eptor,flc_eptor,avg_temp,flc_temp,
     x                 avg_pecek,flc_pecek
               ELSE
                  WRITE(kprint,22050) avg_ptote,flc_ptote,
     x                 avg_pepot,flc_pepot,avg_pecoul,flc_pecoul,
     x                 avg_pnbd,flc_pnbd,avg_enb14,flc_enb14,
     x                 avg_pebnd,flc_pebnd,avg_ebond,flc_ebond,
     x                 avg_ebend,flc_ebend,avg_eitor,flc_eitor,
     x                 avg_eptor,flc_eptor,avg_temp,flc_temp,
     x                 avg_pecek,flc_pecek,avg_temph,flc_temph
               END IF
            END IF
         END IF
         IF(slv_exist .AND. (.NOT. slt_exist) ) THEN
            sum_econf=sum_econf+econf
            sum_ecoul=sum_ecoul+ecoul
            sum_etotpot=sum_etotpot+econf+eslvintra
            sum_eslvint=sum_eslvint+eslvintra
            sum_eebond=sum_eebond+eebond
            sum_eebend=sum_eebend+eebend
            sum_eeptors=sum_eeptors+eeptors
            sum_eeitors=sum_eeitors+eeitors
            sum_tote=sum_tote+tote
            sum_ucek=sum_ucek+ecek
            sum_temp=sum_temp+temp
            sum_tempt=sum_tempt+tempt
            sum_tempr=sum_tempr+tempr
            sum_temph=sum_temph+temph
            sum_rms=sum_rms+rms
            ssm_econf=ssm_econf+econf**2
            ssm_ecoul=ssm_ecoul+ecoul**2
            ssm_etotpot=ssm_etotpot+(econf+eslvintra)**2
            ssm_eslvint=ssm_eslvint+eslvintra**2
            ssm_eebond=ssm_eebond+eebond**2
            ssm_eebend=ssm_eebend+eebend**2
            ssm_eeptors=ssm_eeptors+eeptors**2
            ssm_eeitors=ssm_eeitors+eeitors**2
            ssm_tote=ssm_tote+tote**2
            ssm_ucek=ssm_ucek+ecek**2
            ssm_temp=ssm_temp+temp**2
            ssm_tempt=ssm_tempt+tempt**2
            ssm_tempr=ssm_tempr+tempr**2
            ssm_temph=ssm_temph+temph**2
            ssm_rms=ssm_rms+rms**2
            IF(MOD(mstep,nprop1).EQ.0)THEN
               fnstep=DFLOAT(mstep)
               avg_econf=sum_econf/fnstep
               avg_ecoul=sum_ecoul/fnstep
               avg_etotpot=sum_etotpot/fnstep
               avg_eslvint=sum_eslvint/fnstep
               avg_eebond=sum_eebond/fnstep
               avg_eebend=sum_eebend/fnstep
               avg_eeptors=sum_eeptors/fnstep
               avg_eeitors=sum_eeitors/fnstep
               avg_tote=sum_tote/fnstep
               avg_ucek=sum_ucek/fnstep
               avg_temp=sum_temp/fnstep
               avg_tempt=sum_tempt/fnstep
               avg_tempr=sum_tempr/fnstep
               avg_temph=sum_temph/fnstep
               avg_rmsp=sum_rms/fnstep
               flc_econf=FLUCT(ssm_econf,fnstep,avg_econf)
               flc_ecoul=FLUCT(ssm_ecoul,fnstep,avg_ecoul)
               flc_etotpot=FLUCT(ssm_etotpot,fnstep,avg_etotpot)
               flc_eslvint=FLUCT(ssm_eslvint,fnstep,avg_eslvint)
               flc_eebond=FLUCT(ssm_eebond,fnstep,avg_eebond)
               flc_eebend=FLUCT(ssm_eebend,fnstep,avg_eebend)
               flc_eeptors=FLUCT(ssm_eeptors,fnstep,avg_eeptors)
               flc_eeitors=FLUCT(ssm_eeitors,fnstep,avg_eeitors)
               flc_tote=FLUCT(ssm_tote,fnstep,avg_tote)
               flc_ucek=FLUCT(ssm_ucek,fnstep,avg_ucek)
               flc_temp=FLUCT(ssm_temp,fnstep,avg_temp)
               flc_tempt=FLUCT(ssm_tempt,fnstep,avg_tempt)
               flc_tempr=FLUCT(ssm_tempr,fnstep,avg_tempr)
               flc_temph=FLUCT(ssm_temph,fnstep,avg_temph)
               flc_rms=FLUCT(ssm_rms,fnstep,avg_rmsp)
               WRITE(kprint,40000) fstep
               WRITE(kprint,40100)
               IF(.NOT. (hoover  .OR. thermos)) THEN
                  WRITE(kprint,22000) avg_tote,flc_tote,avg_etotpot
     &                 ,flc_etotpot,avg_econf,flc_econf,avg_ecoul
     &                 ,flc_ecoul,avg_eslvint,flc_eslvint,avg_eebond
     &                 ,flc_eebond,avg_eebend,flc_eebend,avg_eeptors
     &                 ,flc_eeptors,avg_eeitors,flc_eeitors,avg_tempt
     &                 ,flc_tempt,avg_tempr,flc_tempr,avg_temp ,flc_temp
     &                 ,avg_rmsp  ,flc_rms,avg_ucek ,flc_ucek
               ELSE
                  WRITE(kprint,22010) avg_tote,flc_tote,avg_etotpot
     &                 ,flc_etotpot,avg_econf,flc_econf,avg_ecoul
     &                 ,flc_ecoul,avg_eslvint,flc_eslvint,avg_eebond
     &                 ,flc_eebond,avg_eebend,flc_eebend,avg_eeptors
     &                 ,flc_eeptors,avg_eeitors,flc_eeitors,avg_tempt
     &                 ,flc_tempt,avg_tempr,flc_tempr,avg_temp ,flc_temp
     &                 ,avg_rmsp  ,flc_rms,avg_ucek ,flc_ucek,avg_temph
     &                 ,flc_temph
                  
               END IF
            END IF
         END IF
         IF(slv_exist .AND. slt_exist ) THEN
            sum_econf=sum_econf+econf
            sum_ecoul=sum_ecoul+ecoul
            sum_eslvint=sum_eslvint+eslvintra
            sum_enbnd=sum_enbnd+interm
            sum_eebond=sum_eebond+eebond
            sum_eebend=sum_eebend+eebend
            sum_eeptors=sum_eeptors+eeptors
            sum_eeitors=sum_eeitors+eeitors
            sum_tote=sum_tote+tote
            sum_ucek=sum_ucek+ecek
            sum_temp=sum_temp+temp
            sum_tempt=sum_tempt+tempt
            sum_tempr=sum_tempr+tempr
            sum_temph=sum_temph+temph
            sum_temppr=sum_temppr+temppr
            sum_rms=sum_rms+rms
            sum_pecek=sum_pecek+pecek
            sum_pehyd=sum_pehyd+pehyd
            sum_peconf=sum_peconf+peconf
            sum_pecoul=sum_pecoul+pecoul
            sum_percip=sum_percip+percip
            sum_enb14=sum_enb14+enb14
            sum_ebend=sum_ebend+ebend
            sum_ebond=sum_ebond+ebond
            sum_eitor=sum_eitor+eitor
            sum_eptor=sum_eptor+eptor
            sum_pnbd=sum_pnbd+pnbd
            sum_pebnd=sum_pebnd+pebnd
            sum_pepot=sum_pepot+pepot
            sum_ptote=sum_ptote+ptote
            sum_epconf=sum_epconf+epconf
            sum_epcoul=sum_epcoul+epcoul
            
            ssm_econf=ssm_econf+econf**2
            ssm_ecoul=ssm_ecoul+ecoul**2
            ssm_eslvint=ssm_eslvint+eslvintra**2
            ssm_enbnd=ssm_enbnd+interm**2
            ssm_eebond=sum_eebond+eebond**2
            ssm_eebend=sum_eebend+eebend**2
            ssm_eeptors=sum_eeptors+eeptors**2
            ssm_eeitors=sum_eeitors+eeitors**2
            ssm_tote=ssm_tote+tote**2
            ssm_ucek=ssm_ucek+ecek**2
            ssm_temp=ssm_temp+temp**2
            ssm_tempt=ssm_tempt+tempt**2
            ssm_tempr=ssm_tempr+tempr**2
            ssm_temph=ssm_temph+temph**2
            ssm_temppr=ssm_temppr+temppr**2
            ssm_rms=ssm_rms+rms**2
            ssm_pecek=ssm_pecek+pecek**2
            ssm_pehyd=ssm_pehyd+pehyd**2
            ssm_peconf=ssm_peconf+peconf**2
            ssm_pecoul=ssm_pecoul+pecoul**2
            ssm_percip=ssm_percip+percip**2
            ssm_enb14=ssm_enb14+enb14**2
            ssm_ebend=ssm_ebend+ebend**2
            ssm_ebond=ssm_ebond+ebond**2
            ssm_eitor=ssm_eitor+eitor**2
            ssm_eptor=ssm_eptor+eptor**2
            ssm_pnbd=ssm_pnbd+pnbd**2
            ssm_pebnd=ssm_pebnd+pebnd**2
            ssm_pepot=ssm_pepot+pepot**2
            ssm_ptote=ssm_ptote+ptote**2
            ssm_epconf=ssm_epconf+epconf**2
            ssm_epcoul=ssm_epcoul+epcoul**2
            IF(MOD(mstep,nprop1).EQ.0)THEN
               fnstep=DFLOAT(mstep)
               avg_econf=sum_econf/fnstep
               avg_ecoul=sum_ecoul/fnstep
               avg_enbnd=sum_enbnd/fnstep
               avg_eslvint=sum_eslvint/fnstep
               avg_eebond=sum_eebond/fnstep
               avg_eebend=sum_eebend/fnstep
               avg_eeptors=sum_eeptors/fnstep
               avg_eeitors=sum_eeitors/fnstep
               avg_tote=sum_tote/fnstep
               avg_ucek=sum_ucek/fnstep
               avg_temp=sum_temp/fnstep
               avg_tempt=sum_tempt/fnstep
               avg_tempr=sum_tempr/fnstep
               avg_temph=sum_temph/fnstep
               avg_temppr=sum_temppr/fnstep
               avg_rmsp=sum_rms/fnstep
               avg_pecek=sum_pecek/fnstep 
               avg_pehyd=sum_pehyd/fnstep 
               avg_peconf=sum_peconf/fnstep 
               avg_pecoul=sum_pecoul/fnstep 
               avg_percip=sum_percip/fnstep 
               avg_enb14=sum_enb14/fnstep 
               avg_ebend=sum_ebend/fnstep 
               avg_ebond=sum_ebond/fnstep 
               avg_eitor=sum_eitor/fnstep 
               avg_eptor=sum_eptor/fnstep 
               avg_pnbd=sum_pnbd/fnstep 
               avg_pebnd=sum_pebnd/fnstep 
               avg_pepot=sum_pepot/fnstep 
               avg_ptote=sum_ptote/fnstep 
               avg_epconf=sum_epconf/fnstep 
               avg_epcoul=sum_epcoul/fnstep 
               
               flc_econf=FLUCT(ssm_econf,fnstep,avg_econf)
               flc_ecoul=FLUCT(ssm_ecoul,fnstep,avg_ecoul)
               flc_eslvint=FLUCT(ssm_eslvint,fnstep,avg_eslvint)
               flc_enbnd=FLUCT(ssm_enbnd,fnstep,avg_enbnd)
               flc_eebond=FLUCT(ssm_eebond,fnstep,avg_eebond)
               flc_eebend=FLUCT(ssm_eebend,fnstep,avg_eebend)
               flc_eeptors=FLUCT(ssm_eeptors,fnstep,avg_eeptors)
               flc_eeitors=FLUCT(ssm_eeitors,fnstep,avg_eeitors)
               flc_tote=FLUCT(ssm_tote,fnstep,avg_tote)
               flc_ucek=FLUCT(ssm_ucek,fnstep,avg_ucek)
               flc_temp=FLUCT(ssm_temp,fnstep,avg_temp)
               flc_tempt=FLUCT(ssm_tempt,fnstep,avg_tempt)
               flc_tempr=FLUCT(ssm_tempr,fnstep,avg_tempr)
               flc_temph=FLUCT(ssm_temph,fnstep,avg_temph)
               flc_temppr=FLUCT(ssm_temppr,fnstep,avg_temppr)
               flc_rms=FLUCT(ssm_rms,fnstep,avg_rmsp)
               flc_pecek=FLUCT(ssm_pecek,fnstep ,avg_pecek)
               flc_pehyd=FLUCT(ssm_pehyd,fnstep ,avg_pehyd)
               flc_peconf=FLUCT(ssm_peconf,fnstep ,avg_peconf)
               flc_pecoul=FLUCT(ssm_pecoul,fnstep ,avg_pecoul)
               flc_percip=FLUCT(ssm_percip,fnstep ,avg_percip)
               flc_enb14=FLUCT(ssm_enb14,fnstep ,avg_enb14)
               flc_ebend=FLUCT(ssm_ebend,fnstep ,avg_ebend)
               flc_ebond=FLUCT(ssm_ebond,fnstep ,avg_ebond)
               flc_eitor=FLUCT(ssm_eitor,fnstep ,avg_eitor)
               flc_eptor=FLUCT(ssm_eptor,fnstep ,avg_eptor)
               flc_pnbd=FLUCT(ssm_pnbd,fnstep ,avg_pnbd)
               flc_pebnd=FLUCT(ssm_pebnd,fnstep ,avg_pebnd)
               flc_pepot=FLUCT(ssm_pepot,fnstep ,avg_pepot)
               flc_ptote=FLUCT(ssm_ptote,fnstep ,avg_ptote)
               flc_epconf=FLUCT(ssm_epconf,fnstep ,avg_epconf)
               flc_epcoul=FLUCT(ssm_epcoul,fnstep ,avg_epcoul)
               
               WRITE(kprint,40000) fstep
               WRITE(kprint,40100)
               IF(.NOT. (hoover .OR. thermos)) THEN
                  WRITE(kprint,22002) avg_tote,flc_tote,avg_econf
     &                 ,flc_econf,avg_enbnd,flc_enbnd,avg_ecoul
     &                 ,flc_ecoul,avg_eslvint,flc_eslvint,avg_eebond
     &                 ,flc_eebond,avg_eebend,flc_eebend,avg_eeptors
     &                 ,flc_eeptors,avg_eeitors,flc_eeitors,avg_ptote
     &                 ,flc_ptote,avg_pepot,flc_pepot,avg_pecoul
     &                 ,flc_pecoul,avg_peconf,flc_peconf,avg_pehyd
     &                 ,flc_pehyd,avg_pebnd,flc_pebnd,avg_ebond
     &                 ,flc_ebond,avg_ebend,flc_ebend,avg_eitor
     &                 ,flc_eitor,avg_eptor,flc_eptor,avg_epconf
     &                 ,flc_epconf,avg_epcoul,flc_epcoul,avg_temppr
     &                 ,flc_temppr,avg_tempt,flc_tempt,avg_tempr
     &                 ,flc_tempr,avg_temp,flc_temp,avg_ucek,flc_ucek
     &                 ,avg_pecek,flc_pecek
                  
               ELSE
                  WRITE(kprint,22003) avg_tote,flc_tote,avg_econf
     &                 ,flc_econf,avg_enbnd,flc_enbnd,avg_ecoul
     &                 ,flc_ecoul,avg_eslvint,flc_eslvint,avg_eebond
     &                 ,flc_eebond,avg_eebend,flc_eebend,avg_eeptors
     &                 ,flc_eeptors,avg_eeitors,flc_eeitors,avg_ptote
     &                 ,flc_ptote,avg_pepot,flc_pepot,avg_pecoul
     &                 ,flc_pecoul,avg_peconf,flc_peconf,avg_pehyd
     &                 ,flc_pehyd,avg_pebnd,flc_pebnd,avg_ebond
     &                 ,flc_ebond,avg_ebend,flc_ebend,avg_eitor
     &                 ,flc_eitor,avg_eptor,flc_eptor,avg_epconf
     &                 ,flc_epconf,avg_epcoul,flc_epcoul,avg_temppr
     &                 ,flc_temppr,avg_tempt,flc_tempt,avg_tempr
     &                 ,flc_tempr,avg_temp,flc_temp,avg_ucek,flc_ucek
     &                 ,avg_pecek,flc_pecek,avg_temph,flc_temph
               END IF
            END IF
         END IF
         IF(cpress .OR. lberendsen .OR. pressure) THEN
            sum_presst=sum_presst+presst
            sum_press =sum_press+press
            sum_pressc=sum_pressc+pressc
            sum_pressk=sum_pressk+pressk
            sum_volume=sum_volume+volume
            sum_pv=sum_pv+pv
            DO i=1,3
               DO j=1,3
                  sum_st(i,j)=sum_st(i,j)+prt(i,j)
               END DO
            END DO
            sum_co(1,1)=sum_co(1,1)+a
            sum_co(1,2)=sum_co(1,2)+b
            sum_co(1,3)=sum_co(1,3)+c
            sum_co(2,1)=sum_co(2,1)+alpha
            sum_co(2,2)=sum_co(2,2)+bbeta
            sum_co(2,3)=sum_co(2,3)+gamma
            sum_co(3,1)=0.0D0
            sum_co(3,2)=0.0D0
            sum_co(3,3)=0.0D0
            
            ssm_presst=ssm_presst+presst*presst
            ssm_press =ssm_press+press*press
            ssm_pressc=ssm_pressc+pressc*pressc
            ssm_pressk=ssm_pressk+pressk*pressk
            ssm_volume=ssm_volume+volume*volume
            ssm_pv=ssm_pv+pv*pv
            IF(.NOT. lberendsen) THEN
               sum_temppra=sum_temppra+temppra
               ssm_temppra=ssm_temppra+temppra*temppra
            END IF 
            DO i=1,3
               DO j=1,3
                  ssm_st(i,j)=ssm_st(i,j)+prt(i,j)*prt(i,j)
               END DO
            END DO 
            ssm_co(1,1)=ssm_co(1,1)+a*a
            ssm_co(1,2)=ssm_co(1,2)+b*b
            ssm_co(1,3)=ssm_co(1,3)+c*c
            ssm_co(2,1)=ssm_co(2,1)+alpha*alpha
            ssm_co(2,2)=ssm_co(2,2)+bbeta*bbeta
            ssm_co(2,3)=ssm_co(2,3)+gamma*gamma
            ssm_co(3,1)=0.0D0
            ssm_co(3,2)=0.0D0
            ssm_co(3,3)=0.0D0
            IF(MOD(mstep,nprop1) .EQ. 0) THEN
               avg_presst=sum_presst/fnstep
               avg_press =sum_press/fnstep
               avg_pressc=sum_pressc/fnstep
               avg_pressk=sum_pressk/fnstep
               avg_volume=sum_volume/fnstep
               avg_pv=sum_pv/fnstep
               avg_temppra=sum_temppra/fnstep
               DO i=1,3
                  DO j=1,3
                     avg_co(i,j)=sum_co(i,j)/fnstep
                     avg_st(i,j)=sum_st(i,j)/fnstep
                  END DO
               END DO
               flc_presst=FLUCT(ssm_presst,fnstep,avg_presst)
               flc_press =FLUCT(ssm_press ,fnstep,avg_press )
               flc_pressc=FLUCT(ssm_pressc,fnstep,avg_pressc)
               flc_pressk=FLUCT(ssm_pressk,fnstep,avg_pressk)
               flc_volume=FLUCT(ssm_volume,fnstep,avg_volume)
               flc_pv=FLUCT(ssm_pv,fnstep,avg_pv)
               flc_temppra=FLUCT(ssm_temppra,fnstep,avg_temppra)
               DO i=1,3
                  DO j=1,3
                     flc_co(i,j)=FLUCT(ssm_co(i,j),fnstep,avg_co(i,j))
                     flc_st(i,j)=FLUCT(ssm_st(i,j),fnstep,avg_st(i,j))
                  END DO
               END DO
               WRITE(kprint,22100)
               WRITE(kprint,22200)
     &              avg_presst,flc_presst,avg_press,flc_press,avg_pressk
     &              ,flc_pressk,avg_temppra,flc_temppra,avg_volume
     &              ,flc_volume,avg_pv,flc_pv
               WRITE(kprint,22300)
     x              avg_co(1,1),avg_co(1,2),avg_co(1,3),flc_co(1,1)
     &              ,flc_co(1,2),flc_co(1,3),avg_co(2,1),avg_co(2,2)
     &              ,avg_co(2,3),flc_co(2,1),flc_co(2,2),flc_co(2,3)
     &              ,avg_st(1,1),avg_st(1,2),avg_st(1,3),flc_st(1,1)
     &              ,flc_st(1,2),flc_st(1,3),avg_st(2,1),avg_st(2,2)
     &              ,avg_st(2,3),flc_st(2,1),flc_st(2,2),flc_st(2,3)
     &              ,avg_st(3,1),avg_st(3,2),avg_st(3,3),flc_st(3,1)
     &              ,flc_st(3,2),flc_st(3,3)
            END IF
         ELSE
            sum_volume=sum_volume+volume
         END IF
         IF(MOD(mstep,nprop1) .EQ. 0) THEN
            WRITE(kprint,40200)
         END IF
         
         IF(minimize) THEN
            IF(slt_exist .AND. (.NOT. slv_exist )) THEN
               energy=pepot+pv
            END IF
            IF(slv_exist .AND. (.NOT. slt_exist )) THEN
               energy=econf+pv
            END IF
            IF(slt_exist .AND. slv_exist ) THEN
               energy=pepot+econf+epconf+pv
            END IF
         END IF
      END IF

12000 FORMAT(
     &     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' NonBond  = ',f11.2,' Coloumb  = ',f11.2
     &     ,' Bonded   = ',f11.2,/10x,' Stretch  = ',f11.2
     &     ,' Angle    = ',f11.2,' P-Tors   = ',f11.2,
     &     /10x,' I-Tors   = ',f11.2,' TraTemp  = ',f11.2,' ResTemp  = '
     &     ,f11.2,/10x,' TotTemp  = ',f11.2,' RMS      = ',f11.5/)
12010 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' NonBond  = ',f11.2,' Coloumb  = ',f11.2
     &     ,' Bonded   = ',f11.2,/10x,' Stretch  = ',f11.2
     &     ,' Angle    = ',f11.2,' P-Tors   = ',f11.2,
     &     /10x,' I-Tors   = ',f11.2
     &     ,' TraTemp  = ',f11.2,' ResTemp  = ',f11.2,/10x
     &     ,' TotTemp  = ',f11.2,' THoover  = ',f11.2,' RMS      = ',f11
     &     .5/)
13000 FORMAT(
     &     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' NonBond  = ',f11.2,' Coloumb  = ',f11.2
     &     ,' Bonded   = ',f11.2,/10x,' Stretch  = ',f11.2
     &     ,' Angle    = ',f11.2,' P-Tors   = ',f11.2,
     &     /10x,' I-Tors   = ',f11.2,' TraTemp  = ',f11.2,' ResTemp  = '
     &     ,f11.2,/10x,' TotTemp  = ',f11.2,' VarABMD  = ',f11.2,
     &     ' EnrABMD  = ',f11.2/)
13010 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' NonBond  = ',f11.2,' Coloumb  = ',f11.2
     &     ,' Bonded   = ',f11.2,/10x,' Stretch  = ',f11.2
     &     ,' Angle    = ',f11.2,' P-Tors   = ',f11.2,
     &     /10x,' I-Tors   = ',f11.2
     &     ,' TraTemp  = ',f11.2,' ResTemp  = ',f11.2,/10x
     &     ,' TotTemp  = ',f11.2,' THoover  = ',f11.2,' VarABMD  = '
     &     ,f11.2,/10x,' EnrABMD  = ',f11.2/)
12001 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' NonBond  = ',f11.2
     &     ,' Ener14   = ',f11.2,/10x,' Bonded   = ',f11.2
     &     ,' Stretch  = ',f11.2,' Angle    = ',f11.2,/10x
     &     ,' I-Tors   = ',f11.2,' P-Tors   = ',f11.2,' TotTemp  = ',f11
     &     .3,/10x,' ResTemp  = ',f11.2,' TraTemp  = ',f11.2/)
12071 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' NonBond  = ',f11.2
     &     ,' Ener14   = ',f11.2,/10x,' Bonded   = ',f11.2
     &     ,' Stretch  = ',f11.2,' Angle    = ',f11.2,/10x
     &     ,' I-Tors   = ',f11.2,' P-Tors   = ',f11.2,' TotTemp  = ',f11
     &     .3,/10x,' ResTemp  = ',f11.2,' TraTemp  = ',f11.2
     &     ,' ImplSolv = ',f11.2/)
12004 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' NonBond  = ',f11.2
     &     ,' Ener14   = ',f11.2,/10x,' Bonded   = ',f11.2
     &     ,' Stretch  = ',f11.2,' Angle    = ',f11.2,/10x
     &     ,' I-Tors   = ',f11.2,' P-Tors   = ',f11.2,' TotTemp  = ',f11
     &     .3,/10x,' ResTemp  = ',f11.2,' TraTemp  = ',f11.2
     &     ,' VarABMD  = ',f11.4,/10x,' EnrABMD  = ',f11.2/)
12050 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' NonBond  = ',f11.2
     &     ,' Ener14   = ',f11.2,/10x,' Bonded   = ',f11.2
     &     ,' Stretch  = ',f11.2,' Angle    = ',f11.2,/10x
     &     ,' I-Tors   = ',f11.2,' P-Tors   = ',f11.2,' TotTemp  = ',f11
     &     .3,/10x,' Hoover   = ',f11.2,' ResTemp  = ',f11.2
     &     ,' TraTemp  = ',f11.2/)
12054 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' NonBond  = ',f11.2
     &     ,' Ener14   = ',f11.2,/10x,' Bonded   = ',f11.2
     &     ,' Stretch  = ',f11.2,' Angle    = ',f11.2,/10x
     &     ,' I-Tors   = ',f11.2,' P-Tors   = ',f11.2,' TotTemp  = ',f11
     &     .3,/10x,' Hoover   = ',f11.2,' ResTemp  = ',f11.2
     &     ,' TraTemp  = ',f11.2,/10x,' VarABMD  = ',f11.4
     &     ,' EnrABMD  = ',f11.2/)
12006 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' Recipr   = ',f11.2
     &     ,' NonBond  = ',f11.2,/10x,' Ener14   = ',f11.2
     &     ,' Bonded   = ',f11.2,' Stretch  = ',f11.2,/10x
     &     ,' Angle    = ',f11.2,' I-Tors   = ',f11.2,' P-Tors   = ',f11
     &     .3,/10x,' TotTemp  = ',f11.2,' ResTemp  = ',f11.2
     &     ,' TraTemp  = ',f11.2/)
12051 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' TotPot   = '
     &     ,f11.2,/10x,' Coulomb  = ',f11.2,' Recipr   = ',f11.2
     &     ,' NonBond  = ',f11.2,/10x,' Ener14   = ',f11.2
     &     ,' Bonded   = ',f11.2,' Stretch  = ',f11.2,/10x
     &     ,' Angle    = ',f11.2,' I-Tors   = ',f11.2,' P-Tors   = ',f11
     &     .3,/10x,' TotTemp  = ',f11.2,' Hoover   = ',f11.1
     &     ,' ResTemp  = ',f11.2,/10x,' TraTemp  = ',f11.2/)
12002 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' SlvPot   = '
     &     ,f11.2,/10x,' SlvNonbd = ',f11.2,' SlvCoul  = ',f11.2
     &     ,' SlvBnded = ',f11.2,/10x,' SlvStr   = ',f11.2
     &     ,' SlvBend  = ',f11.2,' SlvPtor  = ',f11.2,/10x
     &     ,' SlvItor  = ',f11.2,' SltTot   = ',f11.2,' SltPot   = ',f11
     &     .3,/10x,' SltCoul  = ',f11.2,' SltL-J   = ',f11.2
     &     ,' SltHyd   = ',f11.2,/10x,' SltBond  = ',f11.2
     &     ,' SltStr   = ',f11.2,' SltBen   = ',f11.2,/10x
     &     ,' SltItor  = ',f11.2,' SltPtor  = ',f11.2,' S-SPot   = ',f11
     &     .3,/10x,' S-SCoul  = ',f11.2,' SltTemp  = ',f11.2
     &     ,' TrasTem  = ',f11.2,/10x,' RestTemp = ',f11.2
     &     ,' TotTemp  = ',f11.2/)
12022 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' SlvPot   = '
     &     ,f11.2,/10x,' SlvNonbd = ',f11.2,' SlvCoul  = ',f11.2
     &     ,' SlvBnded = ',f11.2,/10x,' SlvStr   = ',f11.2
     &     ,' SlvBend  = ',f11.2,' SlvPtor  = ',f11.2,/10x
     &     ,' SlvItor  = ',f11.2,' SltTot   = ',f11.2,' SltPot   = ',f11
     &     .3,/10x,' SltCoul  = ',f11.2,' SltL-J   = ',f11.2
     &     ,' SltHyd   = ',f11.2,/10x,' SltBond  = ',f11.2
     &     ,' SltStr   = ',f11.2,' SltBen   = ',f11.2,/10x
     &     ,' SltItor  = ',f11.2,' SltPtor  = ',f11.2,' S-SPot   = ',f11
     &     .3,/10x,' S-SCoul  = ',f11.2,' SltTemp  = ',f11.2
     &     ,' TrasTem  = ',f11.2,/10x,' RestTemp = ',f11.2
     &     ,' TotTemp  = ',f11.2,' ImplSlv =',f11.2)
12003 FORMAT( 
     x     10x,' Tstep    =',f12.1,' Total    = ',f11.2,' SlvPot   = '
     &     ,f11.2,/10x,' SlvNonbd = ',f11.2,' SlvCoul  = ',f11.2
     &     ,' SlvBnded = ',f11.2,/10x,' SlvStr   = ',f11.2
     &     ,' SlvBend  = ',f11.2,' SlvPtor  = ',f11.2,/10x
     &     ,' SlvItor  = ',f11.2,' SltTot   = ',f11.2,' SltPot   = ',f11
     &     .3,/10x,' SltCoul  = ',f11.2,' SltL-J   = ',f11.2
     &     ,' SltHyd   = ',f11.2,/10x,' SltBond  = ',f11.2
     &     ,' SltStr   = ',f11.2,' SltBen   = ',f11.2,/10x
     &     ,' SltItor  = ',f11.2,' SltPtor  = ',f11.2,' S-SPot   = ',f11
     &     .3,/10x,' S-SCoul  = ',f11.2,' SltTemp  = ',f11.2
     &     ,' TrasTem  = ',f11.2,/10x,' RestTemp = ',f11.2
     &     ,' TotTemp  = ',f11.2,' Hoover   = ',f11.2/)
22000 FORMAT( 
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' TotPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' NonBond  = ',f11.2,'+/-',f9.3,' Coulomb  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Bonded   = ',f11.2,'+/-',f9.3,' Stretch  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Angle    = ',f11.2,'+/-',f9.3,' P-Tors   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' I-Tors   = ',f11.2,'+/-',f9.3,' TraTemp  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' ResTemp  = ',f11.2,'+/-',f9.3,' TotTemp  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' RMS      = ',f11.5,'+/-',f9.3,' Kinetic  = ',f11.2
     &     ,'+/-',f9.3)
22010 FORMAT( 
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' TotPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' NonBond  = ',f11.2,'+/-',f9.3,' Coulomb  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' TotBond  = ',f11.2,'+/-',f9.3,' Stretch  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Angle    = ',f11.2,'+/-',f9.3,' P-Tors   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' I-Tors   = ',f11.2,'+/-',f9.3,' TraTemp  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' ResTemp  = ',f11.2,'+/-',f9.3,' TotTemp  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' RMS      = ',f11.5,'+/-',f9.3,' Kinetic  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' THoover  = ',f11.2,'+/-',f9.3)
22002 FORMAT( 
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' SlvPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvNBond = ',f11.2,'+/-',f9.3,' SlvCoul  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvBond  = ',f11.2,'+/-',f9.3,' SlvStr   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvBen   = ',f11.2,'+/-',f9.3,' SlvPtor  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvItor  = ',f11.2,'+/-',f9.3,' SltTot   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltPot   = ',f11.2,'+/-',f9.3,' SltCoul  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltL-J   = ',f11.2,'+/-',f9.3,' SltHyd   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltBond  = ',f11.2,'+/-',f9.3,' SltStr   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltBen   = ',f11.2,'+/-',f9.3,' SltItor  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltPtor  = ',f11.2,'+/-',f9.3,' S-SPot   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' S-SCoul  = ',f11.2,'+/-',f9.3,' SltTemp  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' TrasTem  = ',f11.2,'+/-',f9.3,' RestTemp = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' TotTemp  = ',f11.2,'+/-',f9.3,' SlvKin   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltKin   = ',f11.2,'+/-',f9.3)

22003 FORMAT(
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' SlvPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvNBond = ',f11.2,'+/-',f9.3,' SlvCoul  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvBond  = ',f11.2,'+/-',f9.3,' SlvStr   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvBen   = ',f11.2,'+/-',f9.3,' SlvPtor  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' SlvItor  = ',f11.2,'+/-',f9.3,' SltTot   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltPot   = ',f11.2,'+/-',f9.3,' SltCoul  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltL-J   = ',f11.2,'+/-',f9.3,' SltHyd   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltBond  = ',f11.2,'+/-',f9.3,' SltStr   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltBen   = ',f11.2,'+/-',f9.3,' SltItor  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltPtor  = ',f11.2,'+/-',f9.3,' S-SPot   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' S-SCoul  = ',f11.2,'+/-',f9.3,' SltTemp  = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' TrasTem  = ',f11.2,'+/-',f9.3,' RestTemp = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' TotTemp  = ',f11.2,'+/-',f9.3,' SlvKin   = ',f11.2,
     &     '+/-',f9.3,
     &     /10x,' SltKin   = ',f11.2,'+/-',f9.3,' Hoover   = ',f11.2,
     &     '+/-',f9.3)
22001 FORMAT( 
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' TotPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Coulomb  = ',f11.2,'+/-',f9.3,' NonBond  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Ener14   = ',f11.2,'+/-',f9.3,' Bonded   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Stretch  = ',f11.2,'+/-',f9.3,' Angle    = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' I-Tors   = ',f11.2,'+/-',f9.3,' P-Tors   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' TotTemp  = ',f11.1,'+/-',f9.3,' Kinetic  = ',f11.2
     &     ,'+/-',f9.3)
22050 FORMAT( 
     &     /10x,' Total    = ',f11.2,'+/-',f9.3,' TotPot   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Coulomb  = ',f11.2,'+/-',f9.3,' NonBond  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Ener14   = ',f11.2,'+/-',f9.3,' Bonded   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Stretch  = ',f11.2,'+/-',f9.3,' Angle    = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' I-Tors   = ',f11.2,'+/-',f9.3,' P-Tors   = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' TotTemp  = ',f11.1,'+/-',f9.3,' Kinetic  = ',f11.2
     &     ,'+/-',f9.3,
     &     /10x,' Hoover   = ',f11.1,'+/-',f9.3)

12100 FORMAT(
     x  10x,'                           ---------------              ',
     x     '            ')
12200 FORMAT(
     x     10x,' TotPre   = ',f11.2,' ConPre   = ',f11.2,' KinPre   = '
     &     ,f11.2,/10x,' TmpPre   = ',f11.2,' Volume   = ',f11.2
     &     ,' PV       = ',f11.4)
12300 FORMAT(
     x  '        .....  cell parameters  ....   ',
     x  4x,'      ....       stress    .....'/
     x  1x,'XYZ',1x,3(f10.4,1x),2x,3(f12.4,1x)/
     x  1x,'ABC',1x,3(f10.4,1x),2x,3(f12.4,1x)/
     x  5x,3(2x,8('.'),1x),2x,3(f12.4,1x))

22100 FORMAT(/
     x  10x,'-------------------------- Fluctuating Box -------------',
     x     '-----------'/)
22200 FORMAT(
     x     10x,' TotPre   = ',f11.2,'+/-',f9.3,' ConPre   = ',f11.2
     &     ,'+/-',f9.3,/10x,' KinPre   = ',f11.2,'+/-',f9.3
     &     ,' TmpPre   = ',f11.2,'+/-',f9.2,/10x,' Volume   = ',f11.2
     &     ,'+/-',f9.3,' PV       = ',f11.2,'+/-',f9.3)
22300 FORMAT(
     x  '            .....  cell parameters  ....  ',
     x  '     .....       + / -     .....    '/
     x  9x,3(f10.4,1x),3x,3(f10.7,1x)/
     x  9x,3(f10.4,1x),3x,3(f10.7,1x)/
     x  '            .....        stress     ....   ',
     x  '     .....       + / -     .....'/
     x  6x,3(f11.4,1x),3x,3(f10.5,1x)/
     x  6x,3(f11.4,1x),3x,3(f10.5,1x)/
     x  6x,3(f11.4,1x),3x,3(f10.5,1x))

40000 FORMAT(//80('=')/'=',78x,'='/
     &'=                      Averages over  ',f9.1,
     &' fs of Simulation               ='/'=',78x,'='/80('=')//)
40100 FORMAT(/26x,'**********************************'/
     x        26x,'*                                *'/
     x        26x,'*  Energies in     _Kjoule/mole_ *'/
     x        26x,'*  Temperatures in _Kelvin_      *'/
     x        26x,'*  Pressures in     MPascal      *'/
     x        26x,'*  Volumes   in     A^3          *'/
     x        26x,'*  Distances in     A            *'/
     x        26x,'*  Stress Tensor in Joule/A**3   *'/
     x        26x,'*                                *'/
     x        26x,'**********************************'//)
40200 FORMAT(//80('=')/80('=')///)
50000 FORMAT(10x,' Tstep =',f9.1)
c80200 FORMAT('T ',f11.3,' Angle ',f11.3,' V_ab ',f12.5,' TotPot ',f11.3)
80200 FORMAT('T ',f11.3,' Angle ',f11.3,' V_ab ',f12.5,' NonBond ',f11.3
     &     ,' Tors ',f11.3)
80300 FORMAT('T ',f11.3,' R_g ',f13.6,' V_ab ',f12.5,' NonBond ',f11.3)

      RETURN
      END
