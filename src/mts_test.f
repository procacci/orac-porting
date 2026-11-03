      SUBROUTINE mts_test(efact,ntmtss,fpx_h,fpx_l,fpx_m,fpx_n,ntmtsp
     &     ,timesec,urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m,ucns_h
     &     ,ucns_l,ucns_m,uslvbon,uslvben,uslvtor,uslvitor
     &     ,conf_bnd_slv_n1,coul_bnd_slv,coul_bnd_slv_n1,self_slv,ucek
     &     ,urcp_h,urcp_l,urcp_m,ucop_h,ucop_l,ucop_m,ucnp_h,ucnp_l
     &     ,ucnp_m,ubond,ubend,uitors,uptors,fsbond,fsbend,fsin14,cnb14
     &     ,unb14,cngrp,ungrp,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt,pucek,urcsp_h,urcsp_l,urcsp_m
     &     ,ucosp_h,ucosp_l,ucosp_m,ucnsp_h,ucnsp_l,ucnsp_m,eer_m,eer_l
     &     ,eer_h,nmol,nato,ntap,nstep,ktest,lfirst,maxstp,cpress,volume
     &     ,pext,ucepr,thermos,uceh,hpot,eqdist,eqdist1,force_cost
     &     ,strbonds,nbondsadded,eqang,eqang1,force_ang,strbends
     &     ,nbendsadded,eqdied,eqdied1,force_died,strtors,nitorsadded
     &     ,strvel,ben_vel,tor_vel,strtime0,strtime1,strtim0r
     &     ,strtim1r,rtime,kplotsteer,bondadded,bendadded,torsadded
     &     ,nstepsteer,lstretch,lbend,litor,nplotsteer,time
     &     ,ltest_times,vi,chrge,steer_com,workcom,bussi,ebussi,itprint
     &     ,ifprint)

****Beta Version: Procacci at CECAM   **********************************
*                                                                      *
*     MTS_TEST prints out at each step istantaneous values of          *
*     Energies and forces of selected atoms for all shell              *
*     contributions at the end of a full r-RESPA propagation step.     *    
*     Output is designed to be easily "awked"                          *
*                                                                      *
*     INPUT: all energies and forces                                   *
*     OUTPUT: none                                                     *
*                                                                      *
*                                                                      *
************************************************************************

c======================= DECLARATIONS ==================================

      use meta
      use rem

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      REAL*8 urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m,ucns_h,ucns_l
     &     ,ucns_m,uslvbon,uslvben,ucek,urcp_h,urcp_l,urcp_m ,ucop_h
     &     ,ucop_l,ucop_m,ucnp_h,ucnp_l,ucnp_m,ubond,ubend,uitors
     &     ,uptors,fsbond,fsbend,fsin14,cnb14,unb14,cngrp,ungrp ,pucek
     &     ,urcsp_h,urcsp_l,urcsp_m,ucosp_h,ucosp_l,ucosp_m ,ucnsp_h
     &     ,ucnsp_l,ucnsp_m,eer_m,eer_l,eer_h,conf_bnd_slt_n1
     &     ,coul_bnd_slt,coul_bnd_slt_n1,conf_bnd_slv_n1,coul_bnd_slv
     &     ,coul_bnd_slv_n1,self_slt,self_slv,uslvtor,uslvitor,volume
     &     ,pext,ucepr,uceh,hpot,eqdist(*),eqdist1(*),force_cost(*)
     &     ,eqang(*),eqang1(*),force_ang(*),eqdied(*),eqdied1(*)
     &     ,force_died(*),strtime0,strtime1,rtime,bondadded(*)
     &     ,bendadded(*),torsadded(*),strvel(*),ben_vel(*),tor_vel(*)
     &     ,time,strtim0r,strtim1r,vi(*),chrge(*),workcom,ebussi

      REAL*8 efact,timesec
      REAL*8 fpx_h(*),fpx_l(*),fpx_m(*),fpx_n(*)
      INTEGER  ntmtss(*),ntmtsp(*),nmol,nato,ntap,nstep,ktest,maxstp
     &     ,nbondsadded,nbendsadded,nitorsadded,kplotsteer
     &     ,nstepsteer,lstretch,lbend,litor,nplotsteer,itprint,ifprint
      LOGICAL lfirst,cpress,thermos
      INTEGER strbonds(*),strbends(*),strtors(*)
      LOGICAL ltest_times,steer_com,bussi

c----------------------- LOCAL VARIABLES -------------------------------

      REAL*8 urcs_h0,urcs_l0,urcs_m0,ucos_h0,ucos_l0,ucos_m0,ucns_h0
     &     ,ucns_l0,ucns_m0,uslvbon0,uslvben0,urcp_h0,urcp_l0
     &     ,urcp_m0,ucop_h0,ucop_l0,ucop_m0,ucnp_h0,ucnp_l0,ucnp_m0
     &     ,ubond0,ubend0,uptors0,uitors0,fsbond0,fsbend0,cnb140
     &     ,unb140,cngrp0,ungrp0,ucosp_h0,radfact,wrkbot,wrkbet
     &     ,wrktot,ucosp_l0,ucosp_m0,ucnsp_h0,ucnsp_l0,ucnsp_m0,tim,pv
     &     ,totwrk1,totwrk2,ekinreal,etotreal,aux1
      REAL*8 Nurcs_h,Nurcs_l,Nurcs_m,Nucos_h,Nucos_l,Nucos_m,Nucns_h
     &     ,Nucns_l,Nucns_m,Nuslvbon,Nuslvben,Nurcp_h,Nurcp_l
     &     ,Nurcp_m,Nucop_h,Nucop_l,Nucop_m,Nucnp_h,Nucnp_l,Nucnp_m
     &     ,Nubond,Nubend,Nuptors,Nuitors,Nfsbond,Nfsbend,Nfsin14,Ncnb14
     &     ,Nunb14,Ncngrp,Nungrp,Nurcsp_h,Nurcsp_l,Nurcsp_m,Nucosp_h
     &     ,Nucosp_l,Nucosp_m,Nucnsp_h,Nucnsp_l,Nucnsp_m,ustot,uptot
     &     ,ektot,upstot,utot,fnav,pottot,utot0,wrkbo(500),wrkbe(500)
     &     ,wrkto(500),lambo(500),lambe(500),lamto(500),diffbo(500)
     &     ,diffbe(500),diffto(500),pi,timesteer,cost,forceg,vg,enervi
      REAL*8 fpx_h0(5),fpx_l0(5),fpx_m0(5),fpx_n0(5)
      REAL*8 Nfpx_h(5),Nfpx_l(5),Nfpx_m(5),Nfpx_n(5)
      INTEGER  nav,k,kp,j,i,naus,iaux,ndeg
      logical lskip,lok_st,linsteer,revlinsteer,steer
      save utot0,timesteer,totwrk1,totwrk2,wrkbo,wrkbe,wrkto
     &     ,lok_st
      real*8 un0, un1, unb, rest_n0,rest_n1, rest_nb, rem_nb, n0_energy
     &     ,n1_energy
      real*8 etot, meta_pot,sigma

*     ---- INITIALIZE AND SAVE BETWEEN CALLS ------------------------
      data diffbo, diffbe, diffto  /500*0.D0, 500*0.D0, 500*0.D0/

C -- diffbo, etc are passed to steer_print and used as the value of last
C     call, so they must be saved when leaving and returning to this
C     routine:
      save diffbo, diffbe, diffto 


*==================== EXECUTABLE STATEMENTS ============================

      tim=timesec*nstep

      radfact=0.0174533d0 

      lskip=.false.

      do j=1,5
         if(ntmtsp(j).lt.0) lskip=.true.
      end do
      do j=1,3
         if(ntmtss(j).lt.0) lskip=.true.
      end do

      if(lskip) goto 1005

c----------------------------------------------------------------------------
c     find first average values 
c----------------------------------------------------------------------------

      IF(lfirst) then
         nav=nav+1
         urcs_h0=urcs_h0 +  urcs_h
         urcs_l0=urcs_l0 +  urcs_l
         urcs_m0=urcs_m0 +  urcs_m
         urcp_h0=urcp_h0 +  urcp_h
         urcp_l0=urcp_l0 +  urcp_l
         urcp_m0=urcp_m0 +  urcp_m
         ucos_h0=ucos_h0 +  ucos_h
         ucos_l0=ucos_l0 +  ucos_l
         ucos_m0=ucos_m0 +  ucos_m
         ucop_h0=ucop_h0 +  ucop_h
         ucop_l0=ucop_l0 +  ucop_l
         ucop_m0=ucop_m0 +  ucop_m
         ucns_h0=ucns_h0 +  ucns_h
         ucns_l0=ucns_l0 +  ucns_l
         ucns_m0=ucns_m0 +  ucns_m
         ucnp_h0=ucnp_h0 +  ucnp_h
         ucnp_l0=ucnp_l0 +  ucnp_l
         ucnp_m0=ucnp_m0 +  ucnp_m
         ucosp_h0=ucosp_h0 + ucosp_h
         ucosp_l0=ucosp_l0 + ucosp_l
         ucosp_m0=ucosp_m0 + ucosp_m
         ucnsp_h0=ucnsp_h0 + ucnsp_h
         ucnsp_l0=ucnsp_l0 + ucnsp_l
         ucnsp_m0=ucnsp_m0 + ucnsp_m
         ubond0=ubond0 +   ubond
         ubend0=ubend0 +   ubend
         uptors0=uptors0 + uptors
         uitors0=uitors0 + uitors
         uslvbon0=uslvbon0 +uslvbon
         uslvben0=uslvben0 +uslvben
         fsbond0=fsbond0 +  fsbond
         fsbend0=fsbend0 +  fsbend
         cnb140=cnb140 +    cnb14
         unb140=unb140 +    unb14
         cngrp0=cngrp0 +    cngrp
         ungrp0=ungrp0 +    ungrp  
         k=0
         DO j=1,5
            IF(ntmtsp(j).ne.0) then
               i=ntmtsp(j) 
               k=k+1
               fpx_n0(k)=fpx_n0(i)+fpx_n(i)
               fpx_m0(k)=fpx_m0(i)+fpx_m(i)
               fpx_l0(k)=fpx_l0(i)+fpx_l(i)
               fpx_h0(k)=fpx_h0(i)+fpx_h(i)
            ENDIF  
         END DO   

c----    averages are done: LFIRST set to .FALSE.

         if(nav.GT.15) THEN
            lfirst=.false.
            fnav=1.d0/float(nav)
         ENDIF

         RETURN
      END IF

c------------------------------------------------------------------------
c     print out total energies
c------------------------------------------------------------------------

1005  ustot=0.d0
      ustot=(ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m+urcs_h+urcs_l
     &     +urcs_m+uslvbon+uslvben+uslvtor+uslvitor+conf_bnd_slv_n1
     &     +coul_bnd_slv+coul_bnd_slv_n1+self_slv)*efact/1000.0D0
      uptot=0.d0
      uptot=(urcp_h+urcp_l+urcp_m+ucop_h+ucop_l+ucop_m+ucnp_h+ucnp_l
     &     +ucnp_m+ubond+ubend+uptors+uitors+conf_bnd_slt_n1
     &     +coul_bnd_slt+coul_bnd_slt_n1+self_slt)*efact/1000.d0

      ektot = (ucek+pucek)/1000.d0  
      
      IF(cpress) ektot=ektot+ucepr/1000.0D0
      IF(thermos) ektot=ektot+uceh/1000.0D0

      upstot= 0.d0
      upstot = (urcsp_h+urcsp_l+urcsp_m+ucosp_h+ucosp_l+ucosp_m
     &     +ucnsp_h+ucnsp_l+ucnsp_m)*efact/1000.d0

      ekinreal=ektot


c!!!  vi_test
c      enervi = 0.00
c      do i=1, ntap
c        enervi=enervi+ 0.5*chrge(i)*vi(i)
c      enddo   
c      write(6,*) "ener_recip", efact*enervi/1000.d0 ,efact*eer_m/1000.d0
c     &     ,efact*eer_l/1000.d0,efact*eer_h/1000.d0
c!!!  

#ifdef _MPI_

!------------------------------------------ REM test starts here
      
      if(rem_run) then
        
        if( lfirst .and. rem_ts > 1.d-10 ) then
          
          open(unit=888,file="REM_DIAGNOSTIC",form="formatted")
          
          if(cpress) THEN
            write(888,301)
          else
            write(888,302) 
          endif
          
          lfirst = .false.
        endif
        
        if( mod(tim,rem_ts) < 1.d-10 ) then

! first calculate the three potential terms        
          un0 = ( uslvbon + uslvben + uslvitor +
     &         ubond   + ubend   + uitors ) * efact / 1000.d0
          un1 = ( uslvtor + conf_bnd_slv_n1 + coul_bnd_slv_n1 +
     &         uptors  + conf_bnd_slt_n1 + coul_bnd_slt_n1 ) *
     &         efact / 1000.d0
          unb = ( ucns_h + ucos_h + ucns_l + ucos_l + ucns_m + ucos_m +
     &         urcs_h + urcs_l + urcs_m + coul_bnd_slv + self_slv +
     &         urcp_h + urcp_l + urcp_m + ucop_h + ucop_l + ucop_m +
     &         ucnp_h + ucnp_l + ucnp_m + coul_bnd_slt + self_slt +
     &         urcsp_h + urcsp_l + urcsp_m + ucosp_h + ucosp_l +
     &         ucosp_m + ucnsp_h + ucnsp_l + ucnsp_m + eer_m +
     &         eer_l + eer_h ) * efact / 1000.d0

          if(rem_segment) then  ! calculate the difference between the total potential and the one from the subset
            n0_energy = rem_n0_energy * efact / 1000.d0
            n1_energy = rem_n1_energy * efact / 1000.d0
            rem_nb = ( rem_hnb_energy + rem_lnb_energy +
     &           rem_mnb_energy ) * efact / 1000.d0
            rest_n0 = un0 - n0_energy
            rest_n1 = un1 - n1_energy
            rest_nb = unb - rem_nb
          else
            rest_n0 = 0.d0
            rest_n1 = 0.d0
            rest_nb = 0.d0
          endif

          if(rem_segment) then  ! only the potential terms involving the selected atoms
            rem_pot(1) = n0_energy
            rem_pot(2) = n1_energy
            rem_pot(3) = rem_nb
!            write(6,'(A20,2G15.5)') "mts_test rem_pot N0", n0_energy,
!     &           rem_pot(1)
!            write(6,'(A20,2G15.5)') "mts_test rem_pot N1", n1_energy,
!     &           rem_pot(2)            
          else                  ! use all the potential
            rem_pot(1) = un0
            rem_pot(2) = un1
            rem_pot(3) = unb
          endif

          if(cpress) then

            etot = ( ucek + pucek + ucepr + uceh +
     &           ( pext * volume + hpot ) * efact ) / 1000.d0 +
     &           rem_pot(1) * rem_factor(1) + rest_n0 +
     &           rem_pot(2) * rem_factor(2) + rest_n1 +
     &           rem_pot(3) * rem_factor(3) + rest_nb

            if(mod(timesec*nstep,rem_printd).lt.1.d-10) 
     &           write(888,'(f12.3,6x,i4,1x,10(f12.3,1x))') 
     &           tim,
     &           para_index,
     &           etot,          ! total energy
     &           (ucek + pucek) / 1000.d0, ! kinetic energy
     &           un0 + un1 + unb, ! original (unscaled) potential energy
     &           uceh / 1000.d0, ! thermostat kinetic energy
     &           hpot * efact / 1000.d0, ! thermostat potential energy
     &           ucepr / 1000.0D0, ! barostath kinetic energy
     &           (pext * volume * efact) / 1000.d0 ! barostath potential energy

          else

            etot = (ucek + pucek + uceh + hpot*efact) / 1000.d0 +
     &           rem_pot(1) * rem_factor(1) + rest_n0 +
     &           rem_pot(2) * rem_factor(2) + rest_n1 +
     &           rem_pot(3) * rem_factor(3) + rest_nb

            if(mod(timesec*nstep,rem_printd).lt.1.d-10) 
     &           write(888,'(f12.3,6x,i4,1x,7(f12.3,1x))') 
     &           tim,
     &           para_index, 
     &           etot,          ! total energy
     &           (ucek + pucek) / 1000.d0, ! kinetic energy
     &           un0 + un1 + unb, ! original (unscaled) potential energy
     &           uceh / 1000.d0, ! thermostat kinetic energy
     &           hpot * efact / 1000.d0  ! thermostat potential energy
          endif
        endif

      endif

301   format("#", 3x, "Time    ", 1x, "Replica  ", 9x, "Etot",
     &           9x, "Ekin", 9x, "Epot", 8x, "Ekinh", 8x, "Epoth", 8x,
     &           "EkinP", 11x, "PV")
302   format("#", 3x, "Time    ", 1x, "Replica  ", 9x, "Etot",
     &           9x, "Ekin", 9x, "Epot", 8x, "Ekinh", 8x, "Epoth")
      
!------------------------------------------ REM test ends here
#endif
      
      utot= ektot+ustot+uptot+upstot+(eer_m+eer_l+eer_h)*efact/1000.d0
      pv=0.d0
      IF(cpress) THEN
        pv=(pext*volume*efact)/1000.0D0
        utot=utot+pv
      END IF
      IF(thermos) utot=utot+hpot*efact/1000.0D0
      IF(bussi) utot=utot+ebussi/1000.d0
      IF(meta_run) then 
        utot = utot + meta_pot(meta_rc) * efact / 1000.d0
        utot = utot - (meta_n * meta_height) * efact / 1000.d0
      endif

      pottot = ustot+uptot+upstot
      etotreal=pottot+ekinreal
 
      if(lfirst.and.ltest_times) then 
        write(ktest,290) 
290     FORMAT("#",7x,"time(fs)",6x,"ETOT",3x,6x,"ESLV",3x,6x,"ESLT",3x
     &       ,6x,"ESLV-SLT",6x,"EKIN",3x,6x,"EPTOT",2x,6x,"EREAL",2x,8x
     &       ,"PV",3x,6x,"HPOT",3x,6x,"KINH",3x)
        lfirst=.false.
      END IF
      if(ltest_times) THEN 
        if(mod(nstep,itprint).eq.0) THEN 
          write(ktest,300) tim,utot,ustot,uptot,upstot
     &         ,ekinreal,pottot,etotreal,pv,hpot*efact/1000.d0,uceh/1000
     &         .d0
        END IF
300     FORMAT("time",f12.2,15f13.3)
        if(nstepsteer.eq.0.and.mod(nstep,ifprint).EQ.0) THEN
          if(nbondsadded.gt.0) then
            do i=1,nbondsadded 
              forceg = -2.*force_cost(i)*(bondadded(i)-eqdist(i))
              vg = force_cost(i)*(bondadded(i)-eqdist(i))**2
              write(ktest,450) i,tim,forceg,vg
     &             ,force_cost(i),eqdist(i),bondadded(i)
450           FORMAT('#Bond ',i3,f16.3,5f15.6)
            end do
          endif
          if(nbendsadded.gt.0) then
            do i=1,nbendsadded  
              forceg = -2*force_ang(i)*(bendadded(i)-eqang(i)*radfact)
              vg = force_ang(i)*(bendadded(i)-eqang(i)*radfact)*2
              write(ktest,451) i,tim,forceg,vg
     &             ,force_ang(i),eqang(i),bendadded(i)/radfact
451           FORMAT('#Bend ',i3,f16.3,5f15.6)
            end do
          endif
          if(nitorsadded.gt.0) then
            do i=1,nitorsadded  
              forceg = -2*force_died(i)*(torsadded(i)-eqdied(i)*radfact)
              vg = force_died(i)*(torsadded(i)-eqdied(i)*radfact)**2
              write(ktest,452) i,tim,forceg,vg
     &             ,force_died(i),eqdied(i),torsadded(i)/radfact
452           FORMAT('#Tors ',i3,f16.3,5f15.6)
            end do
          endif
        END IF
      END IF

c------------------------------------------------------------------------
c     print (if any) STEER DYNAMICS data (Jarzynski irrev. wrk function) 
c------------------------------------------------------------------------

1053  if(nstepsteer.eq.0)  lok_st= .false. 

      if(nstepsteer.eq.0.and.((rtime.ge.strtime0).and.(rtime.le
     &     .strtime1))) THEN 
c--     if nplotsteer is not set then do nothing. 
        if(nplotsteer.eq.0) go to 3456 
        utot0=utot
        totwrk1=0.d0
        totwrk2=0.d0
        timesteer=time
        do i=1,nbondsadded
          if(strbonds(i).GT.0) THEN
            lok_st=.true.
            wrkbo(i)=0.d0
            lambo(i)=0.d0
            diffbo(i)=0.d0
          ENDIF 
        end do
        do i=1,nbendsadded
          IF(strbends(i).GT.0) THEN 
            lok_st=.true.
            wrkbe(i)=0.d0
            lambe(i)=0.d0
            diffbe(i)=0.d0
          END IF
        end do
        do i=1,nitorsadded
          if(strtors(i).GT.0) THEN 
            lok_st=.true.
            wrkto(i)=0.d0
            lamto(i)=0.d0
            diffto(i)=0.d0
          endif
        end do
        if(lok_st) THEN 
          write(kplotsteer,1011) 
1011  FORMAT( //"# **** DATA FROM STEERED DYNAMICS *****", //
     &        ,'# The data in this file must be post-processed; '
     &        ,'Lines can be of the following type:  '/
     &        ,' '/
     &        ,'#  a) bond n. I  TIME  R0(T) R  WORK(I) '/
     &        ,'# where: I is the i-the steered bond; TIME is the '
     &        ,' actual time of the simulation; R0(T) is the driven'/
     &        ,'# equilibrium coordinate; R is the actual bond coord' 
     &        ,'inate and WORK(I) is the work done on this coordinate'/ 
     &        ,' '/
     &        ,'#  b) bend n. I  TIME  ALPHA0(T) ALPHA  WORK(I) '/
     &        ,'# where: I is the i-the steered bend; TIME is the '
     &        ,'actual time of the simulation; ALPHA0(T) is the driven'/
     &        ,'# equilibrium coordinate; ALPHA is the actual bend ' 
     &        ,'coordinate and WORK(I) is the work done on this '
     &        ,'coordinate'/ 
     &        ,' '/
     &        ,'#  c) tors n. I  time  ALPHA0(T) ALPHA  WORK(I) '/
     &        ,'# where: I is the i-the steered tors; TIME is the '
     &        ,'actual time of the simulation; ALPHA0(T) is the driven'/
     &        ,'# equilibrium coordinate; ALPHA is the actual bend cord' 
     &        ,'inate and work(i) is the work done on this coordinate'/ 
     &        ,' '/
     &        ,'#  d) WORKfunction T WBO WBE WTO TOTWORK TOTWORK2'/
     &        ,'# where T is the actual time of the simulation ;'
     &        ,' WBO is the total work done on all the stretchings;' /
     &        ,'# WBE is the total work done on all the torsions;' 
     &        ,' TOTWORK is the total work done on the system during'/
     &        ,'# the nonequilibrium process. TOTWORK2 is the total'
     &        ,' work calculated as (total) energy difference. '/
     &        ,'# TOTWORK1 and TOTWORK2 should be equal: a discrepancy'
     &        ,' is due to a bad time integration' //)
        end if
3456    CONTINUE
      END IF

      linsteer= ((rtime.ge.strtime0).and.(rtime.le.strtime1))
      revlinsteer=((rtime.ge.strtim1r).and.(rtime.le.strtim0r))
c     during the reverse run velocity must change sign; work
c     must change sign. Hence we use the same velocity as in the
c     forward path keep for getting -W along the reverse path 
c     The final work is equal (for many realizations) to (about) twice
c     the free energy for the forward path

      steer=linsteer.or.revlinsteer
      if(lok_st.and.steer) THEN 
        nstepsteer=nstepsteer+1 
      END IF
c     do bonds
      if(nplotsteer.gt.0) naus=mod(nstepsteer,nplotsteer)  
c--   Do nothing if nplotsteer is not set
      if(nplotsteer.le.0) go to 3457  
      cost = 1.d0
      aux1=-10.d0
#ifdef _MPI_
      aux1 = rem_ts
      iaux = para_index
#endif
      call steer_print(nbondsadded,force_cost,strvel
     &     ,strbonds,bondadded,eqdist,eqdist1,strtime0,strtime1
     &     ,rtime,wrkbo,lambo,diffbo,timesteer,cost,strtim0r
     &     ,strtim1r)
      wrkbot=0.d0
      if((steer.and.(nbondsadded.gt.0).and.(naus.eq.0)))THEN 
        do  i=1,nbondsadded
          if(strbonds(i).GT.0) THEN 
            wrkbot=wrkbot+ wrkbo(i)
            if(aux1.gt.0.d0) THEN 
              write(kplotsteer,1016) i,rtime,lambo(i)
     &             ,bondadded(i),wrkbo(i),iaux
            ELSE
              write(kplotsteer,1015) i,rtime,lambo(i)
     &             ,bondadded(i),wrkbo(i)
            END IF
1015        FORMAT(" bond n. ", i3,g13.3,4f15.5) 
1016        FORMAT(" bond n. ", i3,g13.3,3f15.5, "  Replica -> ",i4) 
          END IF
        end do
      END IF

c     do bends
      wrkbet=0.d0
      call steer_print(nbendsadded,force_ang,ben_vel
     &     ,strbends,bendadded,eqang,eqang1,strtime0,strtime1,rtime
     &     ,wrkbe,lambe,diffbe,timesteer,radfact,strtim0r,strtim1r)
      if((steer.and.(nbendsadded.gt.0).and.(naus.eq.0)))THEN 

        do i=1,nbendsadded
          if(strbends(i).GT.0) THEN 
            wrkbet=wrkbet+ wrkbe(i)
            if(aux1.gt.0.d0) THEN 
              write(kplotsteer,1026)i,rtime,57.2958
     &             *lambe(i),57.2958*bendadded(i),wrkbe(i),iaux
            ELSE
              write(kplotsteer,1025)i,rtime,57.2958
     &             *lambe(i),57.2958*bendadded(i),wrkbe(i)
            END IF
1025        FORMAT(" bend n. ", i3,g13.2,3f15.5) 
1026        FORMAT(" bend n. ", i3,g13.2,3f15.5, "  Replica -> ",i4)  
          END IF
        end do
      END IF

c     do tors
      wrktot=0.d0
      call steer_print(nitorsadded,force_died,tor_vel
     &     ,strtors,torsadded,eqdied,eqdied1,strtime0,strtime1
     &     ,rtime,wrkto,lamto,diffto,timesteer,radfact,strtim0r
     &     ,strtim1r)
      if(steer.and.(nitorsadded.gt.0).and.(naus.eq.0)) THEN 
        do i=1,nitorsadded
          if(strtors(i).GT.0) THEN 
            wrktot=wrktot+ wrkto(i)
            if(aux1.gt.0.d0) THEN 
              write(kplotsteer,1036)i,rtime,57.2958
     &             *lamto(i),57.2958*torsadded(i),wrkto(i),iaux
            ELSE
              write(kplotsteer,1035)i,rtime,57.2958
     &             *lamto(i),57.2958*torsadded(i),wrkto(i)
            END IF
1035        FORMAT(" tors n. ", i3,g13.2,3f15.5) 
1036        FORMAT(" tors n. ", i3,g13.2,3f15.5, "  Replica -> ",i4)  
          END IF
        end do
      END IF
      
!     do COM 
      if(steer_com) THEN 
        wrktot=wrktot+workcom
      END IF  
      
      if((steer.and.naus.eq.0)) THEN 
        Totwrk1 = wrkbot+wrkbet+wrktot
        if(linsteer) THEN 
          Totwrk2=  Totwrk2 + utot-utot0
        else if(revlinsteer) then
          Totwrk2=  Totwrk2 - (utot-utot0)
        end if
        write(kplotsteer,1045) rtime,wrkbot,wrkbet
     &       ,wrktot,totwrk1,totwrk2
        utot0=utot
      END IF
1045  FORMAT (" WORKfunction",f12.1,5f15.6)   
3457  CONTINUE
      if(lskip) return
      
c------------------------------------------------------------------------
c     shift energies
c------------------------------------------------------------------------

c---  shift energies

      Nurcs_h   = urcs_h   -  urcs_h0 *fnav
      Nurcs_l   = urcs_l   -  urcs_l0 *fnav
      Nurcs_m   = urcs_m   -  urcs_m0 *fnav
      Nurcp_h   = urcp_h   -  urcp_h0 *fnav
      Nurcp_l   = urcp_l   -  urcp_l0 *fnav
      Nurcp_m   = urcp_m   -  urcp_m0 *fnav
      Nucos_h   = ucos_h   -  ucos_h0 *fnav
      Nucos_l   = ucos_l   -  ucos_l0 *fnav
      Nucos_m   = ucos_m   -  ucos_m0 *fnav
      Nucop_h   = ucop_h   -  ucop_h0 *fnav
      Nucop_l   = ucop_l   -  ucop_l0 *fnav
      Nucop_m   = ucop_m   -  ucop_m0 *fnav
      Nucns_h   = ucns_h   -  ucns_h0 *fnav
      Nucns_l   = ucns_l   -  ucns_l0 *fnav
      Nucns_m   = ucns_m   -  ucns_m0 *fnav
      Nucnp_h   = ucnp_h   -  ucnp_h0 *fnav
      Nucnp_l   = ucnp_l   -  ucnp_l0 *fnav
      Nucnp_m   = ucnp_m   -  ucnp_m0 *fnav
      Nucosp_h  = ucosp_h  -  ucosp_h0*fnav
      Nucosp_l  = ucosp_l  -  ucosp_l0*fnav
      Nucosp_m  = ucosp_m  -  ucosp_m0*fnav
      Nucnsp_h  = ucnsp_h -   ucnsp_h0*fnav
      Nucnsp_l  = ucnsp_l -   ucnsp_l0*fnav
      Nucnsp_m  = ucnsp_m -   ucnsp_m0*fnav
      Nubond    = ubond   -   ubond0  *fnav
      Nubend    = ubend   -   ubend0  *fnav
      Nuptors   = uptors -    uptors0 *fnav
      Nuitors   = uitors -    uitors0 *fnav 
      Nuslvbon  = uslvbon -   uslvbon0*fnav
      Nuslvben  = uslvben -   uslvben0*fnav
      Nfsbond   = fsbond  -   fsbond0 *fnav
      Nfsbend   = fsbend  -   fsbend0 *fnav
      Ncnb14    = cnb14   -   cnb140  *fnav
      Nunb14    = unb14   -   unb140  *fnav
      Ncngrp    = cngrp   -   cngrp0  *fnav
      Nungrp    = ungrp   -   ungrp0  *fnav   

c---  shift forces

      kp=0
      DO j=1,5
         IF(ntmtsp(j).ne.0) then
            i=ntmtsp(j) 
            kp=kp+1
            Nfpx_n(kp)=fpx_n(i)-fpx_n0(kp)*fnav
            Nfpx_m(kp)=fpx_m(i)-fpx_m0(kp)*fnav
            Nfpx_l(kp)=fpx_l(i)-fpx_l0(kp)*fnav
            Nfpx_h(kp)=fpx_h(i)-fpx_h0(kp)*fnav
         ENDIF  
      END DO   

c------------------------------------------------------------------------
c     print out "shifted" shells contributions for forces and energies)
c------------------------------------------------------------------------

      DO i=1,kp
        write(ktest,200) i,tim,Nfpx_n(i),Nfpx_m(i)
     &       ,Nfpx_l(i),Nfpx_h(i) 
200      FORMAT('P:Force-',i1,f10.3,4e15.5)
      END DO   
      write(ktest,210) tim,Nurcp_h,Nucop_h,Nucnp_h
      write(ktest,220) tim,Nurcp_l,Nucop_l,Nucnp_l
      write(ktest,230) tim,Nurcp_m,Nucop_m,Nucnp_m
     &     ,Nfsbond,Nfsbend,Nfsin14,cngrp,ungrp,unb14,cnb14
      write(ktest,240) tim,ubond,ubend
210   FORMAT(' P:h-',f10.3,4E15.5)
220   FORMAT(' P:l-',f10.3,4E15.5)
230   FORMAT(' P:m-',f10.3,10E15.5)
240   FORMAT(' P:n-',f10.3,4E15.5)
      write(ktest,110) tim,Nurcs_h,Nucos_h,Nucns_h
      write(ktest,120) tim,Nurcs_l,Nucos_l,Nucns_l
      write(ktest,120) tim,Nurcs_m,Nucos_m,Nucns_m
      write(ktest,140) tim,Nuslvben,Nuslvbon
110   FORMAT(' S:h-',f10.3,4E15.5)
120   FORMAT(' S:l-',f10.3,4E15.5)
130   FORMAT(' S:m-',f10.3,4E15.5)
140   FORMAT(' S:n-',f10.3,4E15.5)
      write(ktest,310) tim,Nurcsp_h,Nucosp_h,Nucnsp_h
      write(ktest,320) tim,Nurcsp_l,Nucosp_l,Nucnsp_l
      write(ktest,330) tim,Nurcsp_m,Nucosp_m,Nucnsp_m
310   FORMAT(' PS:h-',f10.3,4E15.5)
320   FORMAT(' PS:l-',f10.3,4E15.5)
330   FORMAT(' PS:m-',f10.3,4E15.5)
      RETURN
      END
